import json
import logging
from pathlib import Path
import time

# from typing import Generator, Iterator
from matplotlib import pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder
import seaborn as sns
import numpy as np
import xgboost as xgb
from datetime import datetime
from utils import serialize, deserialize, format_duration

# from tqdm.auto import tqdm

from dataset import ByRankGenerator
from taxdb import TaxDB
from database import Database
from dmatstore import DMatStore

LOG = logging.getLogger(__name__)

DEFAULT_PARAMETERS = {
    "max_bin": 256,  # for gpu_hist
    "grow_policy": "depthwise",  # better with big dataset
    "objective": "multi:softmax",  # classification
    "eval_metric": "mlogloss",
    "booster": "gbtree",
    "learning_rate": 0.05,  # eta
    "max_depth": 7,  # 256 columns -> 6 -10
    "min_child_weight": 5,  # prevents creation of unnecessary leaves
    "gamma": 0.1,  # Prunes non-informative branches
    "subsample": 0.8,  # Prevents overfitting by sampling 80% of the data
    "colsample_bytree": 0.8,  # Selects 80% of the features for each tree
    "lambda": 1.0,  # L2 (Ridge) regularization to avoid extreme weight values
    "alpha": 0.5,  # L1 (Lasso) regularization to encourage sparsity
}

EVAL_STORAGE = "eval"
TEST_STORAGE = "test"


class XGBoostModel:
    def __init__(
        self,
        rank: str,
        database: Database,
        # save_path: str | Path,
        batch_size: int,
        # test_batch_count: int,
        # train_batch_count: int | None = None,  # None = max
        # eval_batch_count: int = 0,
        workspace_path: str | Path,
        params: dict | None = None,
        normalize: str | None = None,
        seed: int = 2025,
        taxdb: TaxDB | None = None,
        generator_threads: int = 10,
        sample_balance_factor: float = 0.0,
        batch_balance_factor: float = 0.0,
        min_samples_by_class: int = 1,
        max_buffer_total_size: int | None = None,
    ):
        LOG.debug(f"XGBoostModel({locals()})")
        self._params = params if params is not None else DEFAULT_PARAMETERS
        # self._train_batch_count = train_batch_count
        # self._test_batch_count = test_batch_count
        # self._eval_batch_count = eval_batch_count
        self._seed = seed
        self._rank = rank
        self._normalize = normalize
        self._batch_size = batch_size
        self._taxdb = taxdb
        self._database = database
        self._generator_threads = generator_threads
        self._sample_balance_factor = sample_balance_factor
        self._batch_balance_factor = batch_balance_factor
        self._min_samples_by_class = min_samples_by_class

        self._workspace_path = Path(workspace_path).resolve()
        self._model = None
        self._gen = None
        self._last_batch_id = 0
        self._max_batch_count = None
        self._report = None  # full report
        self._batch_reports = []  # patience
        self._dmat_store = DMatStore(self._workspace_path / "dmatstore")

        # sample generator
        self._gen = ByRankGenerator(
            database=self._database,
            taxdb=self._taxdb,
            rank=self._rank,
            batch_size=self._batch_size,
            normalize=self._normalize,
            seed=self._seed,
            buffer_threads=self._generator_threads,
            sample_balance_factor=self._sample_balance_factor,
            batch_balance_factor=self._batch_balance_factor,
            min_samples_by_class=min_samples_by_class,
            max_buffer_total_size=max_buffer_total_size,
        )
        self._gen.start()

        # labels
        self._labels = self._gen.labels(self._rank, min_samples=min_samples_by_class)
        self._label_encoder = LabelEncoder()
        self._label_encoder.fit(self._labels)

        self._sn_map = {
            tax_id: self._taxdb[tax_id]["ScientificName"] for tax_id in self._labels
        }

        # batches
        self._max_batch_count = self._gen.estimated_batches_count()

        LOG.info(f"XGBoostModel max batches: {self._max_batch_count}")

    def train(
        self,
        save_path: str | Path,
        train_batch_count: int | None = None,
        test_batch_count: int = 1,
        eval_batch_count: int = 1,
        eval_patience: int = 1,
        num_boost_round: int = 100,
    ):
        save_path = Path(save_path).resolve()

        splits = self._get_simple_batch_splits(
            train_batch_count=train_batch_count,
            eval_batch_count=eval_batch_count,
            test_batch_count=test_batch_count,
        )
        report = self._train_and_evaluate(
            eval_patience=eval_patience,
            num_boost_round=num_boost_round,
            splits=splits,
            # train_batch_ids=splits["train"],
            # eval_batch_ids=splits["eval"],
            # test_batch_ids=splits["test"],
        )

        self.save(save_path)
        self._generate_report(path=save_path / "report", report=report)

    def kfold(
        self,
        k: int = 5,
        train_batch_count: int | None = None,
        eval_batch_count: int = 1,
        eval_patience: int = 1,
        num_boost_round: int = 100,
    ):
        ksplits = self._get_kfold_batch_splits(
            k, train_batch_count=train_batch_count, eval_batch_count=eval_batch_count
        )
        for splits in ksplits:
            self._train_and_evaluate(
                eval_patience=eval_patience,
                num_boost_round=num_boost_round,
                splits=splits,
            )

    def _train_and_evaluate(
        self, eval_patience: int, num_boost_round: int, splits: dict
    ) -> dict:
        report = {}
        batchs_report = []
        self._model = None
        train_batch_ids = splits["train"]
        eval_batch_ids = splits["eval"]
        test_batch_ids = splits["test"]

        report["start_dt"] = datetime.now()
        # params
        params = self._params.copy()
        if "seed" not in params:
            params["seed"] = self._seed
        params["num_class"] = len(self._labels)

        report["params"] = params
        report["train"] = {"total_duration": 0, "batches": []}
        report["num_boost_round"] = num_boost_round
        report["splits"] = splits
        start_total_time = time.time()

        for i, batch_id in enumerate(train_batch_ids):
            # 1 - TRAIN
            LOG.debug(f"Train : {i+1} / {len(train_batch_ids)} - batch ID: {batch_id}")

            dtrain = self._get_batch(batch_id)
            start_batch_time = time.time()
            batch_report = {
                "batch_duration": 0,
                "get_batch_duration": 0,
                "train_duration": 0,
                "eval_duration": 0,
                "report": None,
                "eval_score": None,
            }

            y = dtrain.get_label().astype(int)
            y_encoded = self._label_encoder.transform(y)
            dtrain_encoded = xgb.DMatrix(dtrain.get_data(), label=y_encoded)
            batch_report["get_batch_duration"] = time.time() - start_batch_time
            start_train_time = time.time()
            self._model = xgb.train(
                params,
                dtrain_encoded,
                num_boost_round=num_boost_round,
                xgb_model=self._model,
            )
            batch_report["train_duration"] = time.time() - start_train_time
            if eval_batch_ids:
                # 2 - EVAL
                start_eval_time = time.time()
                eval_report = self._evaluate(batch_ids=eval_batch_ids)
                score = eval_report["score"]
                LOG.debug(f"Evaluation score: {score:.3f}")
                self.save(self._workspace_path / "batches " / str(len(batchs_report)))
                prev_scores = [br["score"] for br in batchs_report]
                prev_scores_str = ", ".join(f"{ps:.2f}" for ps in prev_scores)
                LOG.debug(
                    f"Batch score: {score:.2f}, previous scores: {prev_scores_str}"
                )
                batchs_report.append(eval_report)

                batch_report.update(
                    {
                        "eval_duration": time.time() - start_eval_time,
                        "report": eval_report,
                        "eval_score": score,
                    }
                )

                if len(batchs_report) >= eval_patience:
                    patience_prev_scores = prev_scores[-eval_patience:]
                    min_improvement = 0.001  # TODO: conf
                    if score < all(
                        prev_score > score + min_improvement
                        for prev_score in patience_prev_scores
                    ):
                        LOG.info(
                            f"Early stopping: no improvement in the last {eval_patience} evaluations"
                        )
                        break

            report["train"]["batches"].append(batchs_report)

        LOG.debug("Loading best model...")
        prev_scores = [br["score"] for br in batchs_report]
        max_index, max_value = max(enumerate(prev_scores), key=lambda x: x[1])
        LOG.debug(f"Loading best model: BATCH {max_index} => {max_value:.3f}")
        self.load(self._workspace_path / "batches " / str(max_index))

        report["train"]["total_duration"] = time.time() - start_total_time
        report["train"]["batches"] = batchs_report
        LOG.debug("Model trained")

        # 3 - TEST
        start_test_time = time.time()
        report["test"] = self._evaluate(batch_ids=test_batch_ids)
        LOG.debug(f"Test score: {score:.3f}")
        report["test"]["eval_duration"] = time.time() - start_test_time
        report["test"]["best_batch"] = ({"index": max_index, "score": max_value},)
        report["end_dt"] = datetime.now()
        return report

    def _evaluate(self, batch_ids: list[int]) -> dict:
        """Evaluate the model (eval or test)"""

        LOG.debug("Evaluating model...")

        all_y_true_encoded = []
        all_y_pred_encoded = []
        report = {}

        # evaluate
        for i, batch_id in enumerate(batch_ids):
            deval = self._get_batch(batch_id)
            LOG.debug(
                f"Evaluation batch {i + 1} / {len(batch_ids)}) - batch ID: {batch_id}"
            )

            y_true = deval.get_label().astype(int)
            y_true_encoded = self._label_encoder.transform(y_true)
            dteval_encoded = xgb.DMatrix(deval.get_data(), label=y_true_encoded)
            y_pred_encoded = self._model.predict(dteval_encoded).astype(int)

            all_y_true_encoded.extend(y_true_encoded)
            all_y_pred_encoded.extend(y_pred_encoded)

        y_true = self._label_encoder.inverse_transform(all_y_true_encoded)
        y_pred = self._label_encoder.inverse_transform(all_y_pred_encoded)

        # scientific names
        y_true = [self._sn_map[tax_id] for tax_id in y_true]
        y_pred = [self._sn_map[tax_id] for tax_id in y_pred]

        LOG.debug("Model evaluated, getting report...")

        # generate classification report and confusion matrix
        sorted_labels = sorted(self._sn_map.values())
        report["classification_report"] = classification_report(
            y_true,
            y_pred,
            labels=sorted_labels,
            output_dict=True,
        )
        report["confusion_matrix"] = confusion_matrix(
            y_true, y_pred, labels=sorted_labels
        )
        report["score"] = self._score(report)

        LOG.debug("Model evaluated")
        return report

    def _generate_report(self, path: Path, report: dict):
        report_lines = []
        dt_format = "%Y-%m-%d %H:%M:%S"

        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        # 1 - main report
        report_lines.append("=== Training / Evaluation Report ===")
        report_lines.append(f"Rank: {self._rank}")
        report_lines.append(f"Start Date: {report['start_dt'].strftime(dt_format)}")
        report_lines.append(f"End Date: {report['end_dt'].strftime(dt_format)}")
        total_duration_s = (report["end_dt"] - report["start_dt"]).total_seconds()
        report_lines.append(f"Total Duration: {format_duration(total_duration_s)}")
        report_lines.append("")

        report_lines.append("\n=== Parameters ===")
        report_lines.append(f"Database: {self._database.get_db_path()}")
        report_lines.append(f"Boosting rounds: {report['num_boost_round']}")
        report_lines.append("Parameters:")
        report_lines.append(json.dumps(report["params"], indent=4))
        report_lines.append("")

        report_lines.append("\n=== Training Data ===")
        report_lines.append(f"Samples per Batch: {self._batch_size}")
        report_lines.append(f"Sample balance factor: {self._sample_balance_factor}")
        report_lines.append(f"Batch balance factor: {self._batch_balance_factor}")
        report_lines.append(f"Min samples per class: {self._min_samples_by_class}")
        report_lines.append(f"Available Batches: {self._max_batch_count}")
        report_lines.append(f"Actually used Batches: {self._last_batch_id}")
        report_lines.append(
            f"Training Batches: {len(report['train']['batches'])} (max: {len(report['splits']['train'])})"
        )
        report_lines.append(f"Evaluation Batches: {len(report['splits']['eval'])}")
        report_lines.append(f"Test Batches: {len(report['splits']['test'])}")
        report_lines.append(f"Normalization: {self._normalize}")
        report_lines.append(f"Seed: {self._seed}")
        report_lines.append("")

        report_lines.append("\n=== Classification Report ===")
        for label, metrics in report["test"]["classification_report"].items():
            if isinstance(metrics, dict):
                report_lines.append(f"Label : {label}")
                for metric_name, value in metrics.items():
                    report_lines.append(f"  {metric_name} : {value:.4f}")
            else:
                report_lines.append(f"{label} : {metrics:.4f}")

        report_lines.append("\n=== Confusion Matrix ===")
        report_lines.append(
            self._confusion_matrix_to_ascii(report["test"]["confusion_matrix"])
        )
        self._plot_and_save_confusion_matrix(
            report["test"]["confusion_matrix"],
            path / "confusion_matrix.png",
        )

        report_str = "\n".join(report_lines)

        report_file = path / "report.txt"
        report_file.write_text(report_str, encoding="utf-8")

        # 2 - tax_id report
        tid_report = self._get_tax_id_report()
        tid_report_lines = []

        non_skipped_ranks = [info for info in tid_report if not info["skipped"]]
        skipped_ranks = [info for info in tid_report if info["skipped"]]

        for info in non_skipped_ranks:
            rank_name = info["name"]
            rank_tid = info["tax_id"]
            tid_report_lines.append(
                f"{self._rank}: {rank_name} [{rank_tid}] (Total count: {info['count']})"
            )
            tid_report_lines.extend(
                f"  - {tax_info['name']} [{tax_id}] (count: {tax_info['count']})"
                for tax_id, tax_info in info["tax_ids"].items()
            )

        if skipped_ranks:
            tid_report_lines.append("\nSkipped Ranks:")
            for info in skipped_ranks:
                rank_name = info["name"]
                rank_tid = info["tax_id"]
                tid_report_lines.append(
                    f"Rank: {rank_name} [{rank_tid}] (Total count: {info['count']})"
                )
                tid_report_lines.extend(
                    f"  - {tax_info['name']} [{tax_id}] (count: {tax_info['count']})"
                    for tax_id, tax_info in info["tax_ids"].items()
                )

        tid_report_str = "\n".join(tid_report_lines)
        tid_report_file = path / "tax_report.txt"
        tid_report_file.write_text(tid_report_str, encoding="utf-8")

        # 3 - raw report (json)
        raw_report_path = path / "raw_report.txt"
        raw_report_path.write_text(
            json.dumps(
                report,
                indent=4,
                default=lambda obj: (
                    obj.isoformat() if isinstance(obj, datetime) else None
                ),
            )
        )

    def _get_tax_id_report(self) -> list[dict]:
        tid_by_rank = self._gen._get_tax_ids_by_rank(self._rank)
        gen_tid_by_rank = self._gen._tax_ids_by_rank
        counts = self._database.counts()
        report = []

        for rank_tid, tax_ids in tid_by_rank.items():
            rank_report = {
                "name": self._taxdb[rank_tid].get("ScientificName"),
                "tax_id": rank_tid,
                "skipped": rank_tid not in gen_tid_by_rank,
            }

            tax_ids_report = {}
            for tax_id in tax_ids:
                tax_ids_report[tax_id] = {
                    "name": self._taxdb[tax_id].get("ScientificName", "Unknown"),
                    "count": counts.get(tax_id, 0),
                    "tax_id": tax_id,
                }

            rank_report["tax_ids"] = tax_ids_report
            rank_report["count"] = sum(
                tax_id["count"] for tax_id in tax_ids_report.values()
            )
            report.append(rank_report)

        return report

    def load(self, path: str | Path) -> None:
        """Load model from file."""
        path = Path(path).resolve()
        LOG.debug(f"Loading model: {path}...")
        model_path = path / "model.bin"
        label_path = path / "labels.pkl"
        label_encoder_path = path / "labels_encoder.pkl"
        params_path = path / "params.json"
        for f in (model_path, label_path, label_encoder_path, params_path):
            if not f.is_file:
                self._gen.stop()
                raise FileNotFoundError(f)

        self._model = xgb.Booster()
        self._model.load_model(str(model_path))
        self._labels = deserialize(label_path)
        self._label_encoder = deserialize(label_encoder_path)
        self._params = json.loads(params_path.read_text())

    def save(self, path: str | Path) -> None:
        """Save model to directory."""
        path = Path(path).resolve()
        LOG.debug(f"Saving model: {path}...")
        path.mkdir(parents=True, exist_ok=True)
        model_path = path / "model.bin"
        label_path = path / "labels.pkl"
        label_encoder_path = path / "labels_encoder.pkl"
        params_path = path / "params.json"

        if self._model is not None:
            self._model.save_model(str(model_path))
            serialize(self._labels, label_path)
            serialize(self._label_encoder, label_encoder_path)
            params_path.write_text(json.dumps(self._params, indent=4))
        else:
            raise ValueError("Model is not trained or loaded yet.")

    def stop(self) -> None:
        self._gen.stop()

    def _score(self, report: dict, lambda_penalty: float = 0.5) -> float:
        """Custom metric: f1 avg, but add penalty according to rank position"""
        f1_macro_avg = report["classification_report"]["macro avg"]["f1-score"]

        f1_scores = [
            v["f1-score"]
            for k, v in report["classification_report"].items()
            if k not in ["accuracy", "macro avg", "weighted avg"]
        ]

        # normalization
        f1_min, f1_max = min(f1_scores), max(f1_scores)
        if f1_max > f1_min:  # Éviter division par zéro
            f1_normalized = [(f - f1_min) / (f1_max - f1_min) for f in f1_scores]
        else:
            f1_normalized = [1] * len(f1_scores)  # Cas où tous les F1 sont identiques

        # apply penalty according to rank position
        penalty = sum((1 - f) for f in f1_normalized) / len(f1_scores)

        return f1_macro_avg - lambda_penalty * penalty

    def _confusion_matrix_to_ascii(self, conf_matrix: list) -> str:
        sorted_labels = sorted(self._sn_map.values())
        max_label_length = max(len(label) for label in sorted_labels)
        header = f"{'':>{max_label_length+2}} " + " ".join(
            f"{label:>{max_label_length}}" for label in sorted_labels
        )
        lines = [header]

        for i, label in enumerate(sorted_labels):
            line = f"{label:{max_label_length}}: " + " ".join(
                f"{conf_matrix[i, j]:{max_label_length}}"
                for j in range(len(sorted_labels))
            )
            lines.append(line)
        return "\n".join(lines)

    def _plot_and_save_confusion_matrix(
        self, conf_matrix: list, file_path: Path
    ) -> None:
        sorted_labels = sorted(self._sn_map.values())
        _, ax = plt.subplots(figsize=(10, 7))
        sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", ax=ax, cbar=False)

        ax.set_xlabel("Predicted Labels")
        ax.set_ylabel("True Labels")
        ax.set_title("Confusion Matrix")

        ax.set_xticks(np.arange(len(sorted_labels)) + 0.5)
        ax.set_yticks(np.arange(len(sorted_labels)) + 0.5)
        ax.set_xticklabels(sorted_labels, rotation=45, ha="right")
        ax.set_yticklabels(sorted_labels, rotation=0)

        plt.tight_layout()
        plt.savefig(file_path)
        plt.close()

    def _get_simple_batch_splits(
        self,
        train_batch_count: int | None,
        test_batch_count: int,
        eval_batch_count: int,
    ) -> dict[str, list[int]]:
        min_required_batches = 1 + eval_batch_count + test_batch_count
        if train_batch_count is None:
            min_required_batches += 1
        else:
            min_required_batches += train_batch_count

        if self._max_batch_count < min_required_batches:
            raise ValueError(
                f"Not enough batches available for the specified split ({self._max_batch_count})"
            )

        test_indices = list(range(1, test_batch_count + 1))
        eval_indices = list(
            range(
                test_batch_count + 1,
                test_batch_count + eval_batch_count + 1,
            )
        )

        train_start_index = test_batch_count + eval_batch_count + 1

        if train_batch_count is None:
            train_indices = list(range(train_start_index, self._max_batch_count + 1))
        else:
            train_indices = list(
                range(train_start_index, train_start_index + train_batch_count)
            )

        return {"train": train_indices, "eval": eval_indices, "test": test_indices}

    def _get_kfold_batch_splits(
        self, k: int, train_batch_count: int | None, eval_batch_count: int
    ):
        if self._max_batch_count < k + eval_batch_count:
            raise ValueError(
                f"Not enough batches available for the specified k-fold split ({self._max_batch_count})"
            )

        splits = []
        all_batches = list(range(1, self._max_batch_count + 1))
        fold_size = (len(all_batches) + k - 1) // k

        for i in range(k):
            start_index = i * fold_size
            end_index = min(start_index + fold_size, len(all_batches))

            test_indices = all_batches[start_index:end_index]
            remaining_indices = [j for j in all_batches if j not in test_indices]
            eval_index = remaining_indices[:eval_batch_count]

            if train_batch_count is None:
                train_indices = remaining_indices[eval_batch_count:]
            else:
                train_indices = remaining_indices[
                    eval_batch_count : eval_batch_count + train_batch_count
                ]

            splits.append(
                {"train": train_indices, "eval": eval_index, "test": test_indices}
            )

        return splits

    def _get_next_batch(self) -> xgb.DMatrix:
        self._last_batch_id += 1
        try:
            return next(self._gen.get())
        except StopIteration:
            LOG.exception(
                f"No more batch available: {self._last_batch_id} / {self._max_batch_count}"
            )
            self._last_batch_id -= 1
            raise

    def _get_batch(self, batch_id: int) -> xgb.DMatrix:
        while not self._dmat_store.has_dmat(batch_id):
            batch = self._get_next_batch()
            self._dmat_store[self._last_batch_id] = batch

        return self._dmat_store[batch_id]
