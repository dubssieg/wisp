import json
import logging
from pathlib import Path
import time
from typing import Generator
from matplotlib import pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder
import seaborn as sns
import numpy as np
import xgboost as xgb
from datetime import datetime
from utils import serialize, deserialize, format_duration, format_size

from dataset import ByRankGenerator
from api import API
from database import Database

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
        save_path: str | Path,
        batch_size: int,
        test_batch_count: int,
        train_batch_count: int | None = None,  # None = max
        eval_batch_count: int = 0,
        params: dict | None = None,
        normalize: str | None = None,
        seed: int = 2025,
        api: API | None = None,
        generator_threads: int = 10,
        sample_balance_factor: float = 0.0,
        batch_balance_factor: float = 0.0,
        min_samples_by_class: int = 1,
        max_buffer_total_size: int | None = None,
    ):
        LOG.debug(f"XGBoostModel({locals()})")
        self._params = params if params is not None else DEFAULT_PARAMETERS
        self._train_batch_count = train_batch_count
        self._test_batch_count = test_batch_count
        self._eval_batch_count = eval_batch_count
        self._seed = seed
        self._rank = rank
        self._normalize = normalize
        self._batch_size = batch_size
        self._api = api
        self._database = database
        self._generator_threads = generator_threads
        self._sample_balance_factor = sample_balance_factor
        self._batch_balance_factor = batch_balance_factor
        self._save_path = save_path = Path(save_path).resolve()
        self._model = None
        self._gen = None
        self._last_batch_id = None
        self._max_batch_count = None
        self._report = None  # full report
        self._batch_reports = []  # patience

        # sample generator
        self._gen = ByRankGenerator(
            database=self._database,
            api=self._api,
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
            tax_id: self._api[tax_id]["ScientificName"] for tax_id in self._labels
        }

        # batches
        self._max_batch_count = self._gen.estimated_batches_count()
        if self._train_batch_count is None:
            self._train_batch_count = max(
                1,
                self._max_batch_count - self._test_batch_count - self._eval_batch_count,
            )

        needed_batches = (
            self._train_batch_count + self._eval_batch_count + self._test_batch_count
        )
        if needed_batches > self._max_batch_count:
            raise ValueError(
                f"not enough batches available: {needed_batches} > {self._max_batch_count}"
            )

        LOG.info(
            f"XGBoostModel batches (size: {self._batch_size}): max: {self._max_batch_count}, train: {self._train_batch_count}, eval: {self._eval_batch_count}, test: {self._test_batch_count}"
        )

    def _store_eval_and_test_batches(self) -> None:
        # store eval batches
        if self._eval_batch_count:
            LOG.debug(
                f"Generating and storing {self._eval_batch_count} eval batches..."
            )
            self._store_batches(
                batch_count=self._eval_batch_count, storage_name=EVAL_STORAGE
            )
            LOG.debug(f"{self._eval_batch_count} eval batches stored")
        LOG.debug(f"Generating and storing {self._test_batch_count} test batches...")
        self._store_batches(
            batch_count=self._test_batch_count, storage_name=TEST_STORAGE
        )
        LOG.debug(f"{self._test_batch_count} test batches stored")

    def train(
        self,
        eval_patience: int = 1,
        num_boost_round: int = 100,
    ) -> None:
        """Train a model"""

        self._store_eval_and_test_batches()

        LOG.debug("Training model...")
        self._report = {}
        self._report["start_dt"] = datetime.now()

        # params
        params = self._params.copy()
        if "seed" not in params:
            params["seed"] = self._seed
        params["num_class"] = len(self._labels)

        self._report["params"] = params
        self._report["train"] = {"total_duration": 0, "batches": []}
        self._report["num_boost_round"] = num_boost_round
        start_total_time = time.time()

        # model
        self._model = None
        for i in range(self._train_batch_count):
            start_batch_time = time.time()
            batch_report = {
                "batch_duration": 0,
                "get_batch_duration": 0,
                "train_duration": 0,
                "eval_duration": 0,
                "report": None,
                "eval_score": None,
            }

            dtrain = self._get_next_batch()

            LOG.debug(f"Train batch {i + 1} / {self._train_batch_count}")

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

            if self._eval_batch_count:
                start_eval_time = time.time()
                report = self.evaluate(storage_name=EVAL_STORAGE)
                score = report["score"]
                LOG.debug(f"Evaluation score: {score:.3f}")
                self.save(f"batch_models/{len(self._batch_reports)}")

                prev_scores = [br["score"] for br in self._batch_reports]
                prev_scores_str = ", ".join(f"{ps:.2f}" for ps in prev_scores)
                LOG.debug(
                    f"Batch score: {score:.2f}, previous scores: {prev_scores_str}"
                )
                self._batch_reports.append(report)

                batch_report.update(
                    {
                        "eval_duration": time.time() - start_eval_time,
                        "report": report,
                        "eval_score": score,
                    }
                )

                if len(self._batch_reports) >= eval_patience:
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

            self._report["train"]["batches"].append(batch_report)

        self._load_best_model()
        self._report["train"]["total_duration"] = time.time() - start_total_time
        LOG.debug("Model trained")

    def evaluate(
        self, batches: list[xgb.DMatrix] | None = None, storage_name: str | None = None
    ) -> dict:
        """Evaluate the model on the next available batches."""

        LOG.debug("Evaluating model...")
        if self._model is None:
            raise ValueError("Model has not been trained yet.")

        all_y_true_encoded = []
        all_y_pred_encoded = []
        report = {}

        if batches:
            if isinstance(batches, xgb.DMatrix):
                batches = [batches]
            batch_iterator = iter(batches)
        else:
            if storage_name is None:
                storage_name = TEST_STORAGE
            batch_iterator = self._stored_batches_generator(storage_name)

        # evaluate
        for i, dtest in enumerate(batch_iterator):
            LOG.debug(f"Evaluation batch {i + 1} / {self._test_batch_count})")

            y_true = dtest.get_label().astype(int)
            y_true_encoded = self._label_encoder.transform(y_true)
            dtest_encoded = xgb.DMatrix(dtest.get_data(), label=y_true_encoded)
            y_pred_encoded = self._model.predict(dtest_encoded).astype(int)

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
        if storage_name == TEST_STORAGE:
            self._report["evaluation"] = report
            self._report["end_dt"] = datetime.now()
        return report

    def generate_report(self):
        report_lines = []
        report = self._report
        dt_format = "%Y-%m-%d %H:%M:%S"

        report_dir = self._save_path / "report"
        report_dir.mkdir(parents=True, exist_ok=True)

        # basic
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

        # training
        report_lines.append("\n=== Training Data ===")
        report_lines.append(f"Samples per Batch: {self._batch_size}")
        report_lines.append(f"Sample balance factor: {self._sample_balance_factor}")
        report_lines.append(f"Batch balance factor: {self._batch_balance_factor}")
        report_lines.append(f"Min samples per class: {self._batch_balance_factor}")
        report_lines.append(f"Available Batches: {self._max_batch_count}")
        report_lines.append(f"Actually used Batches: {self._last_batch_id + 1}")
        train_batches = (
            self._last_batch_id + 1 - self._eval_batch_count - self._test_batch_count
        )
        report_lines.append(
            f"Training Batches: {train_batches} (max: {self._train_batch_count})"
        )

        report_lines.append(f"Evaluation Batches: {self._eval_batch_count}")
        report_lines.append(f"Test Batches: {self._test_batch_count}")
        report_lines.append(f"Normalization: {self._normalize}")
        report_lines.append(f"Seed: {self._seed}")
        report_lines.append("")

        report_lines.append("\n=== Classification Report ===")
        for label, metrics in report["evaluation"]["classification_report"].items():
            if isinstance(metrics, dict):
                report_lines.append(f"Label : {label}")
                for metric_name, value in metrics.items():
                    report_lines.append(f"  {metric_name} : {value:.4f}")
            else:
                report_lines.append(f"{label} : {metrics:.4f}")

        report_lines.append("\n=== Confusion Matrix ===")
        report_lines.append(
            self._confusion_matrix_to_ascii(report["evaluation"]["confusion_matrix"])
        )
        self._plot_and_save_confusion_matrix(
            report["evaluation"]["confusion_matrix"],
            report_dir / "confusion_matrix.png",
        )

        report_str = "\n".join(report_lines)

        report_file = report_dir / "report.txt"
        report_file.write_text(report_str, encoding="utf-8")

    def load(self, sub_dir: str | None = None) -> None:
        """Load model from file."""
        save_path = self._save_path
        if sub_dir:
            save_path /= sub_dir
        LOG.debug(f"Loading model: {save_path}")
        model_path = save_path / "model.bin"
        label_path = save_path / "labels.pkl"
        label_encoder_path = save_path / "labels_encoder.pkl"
        params_path = save_path / "params.json"
        for f in (model_path, label_path, label_encoder_path, params_path):
            if not f.is_file:
                raise FileNotFoundError(f)

        self._model = xgb.Booster()
        self._model.load_model(str(model_path))
        self._labels = deserialize(label_path)
        self._label_encoder = deserialize(label_encoder_path)
        self._params = json.loads(params_path.read_text())

    def _load_best_model(self) -> None:
        LOG.debug("Loading best model...")
        prev_scores = [br["score"] for br in self._batch_reports]
        max_index, max_value = max(enumerate(prev_scores), key=lambda x: x[1])
        LOG.debug(f"Loading best model: BATCH {max_index} => {max_value:.3f}")
        self.load(f"batch_models/{max_index}")

    def save(self, sub_dir: str | None = None) -> None:
        """Save model to directory."""
        save_path = self._save_path
        if sub_dir:
            save_path /= sub_dir
            LOG.debug(f"Saving model: {save_path}")
        save_path.mkdir(parents=True, exist_ok=True)
        model_path = save_path / "model.bin"
        label_path = save_path / "labels.pkl"
        label_encoder_path = save_path / "labels_encoder.pkl"
        params_path = save_path / "params.json"

        if self._model is not None:
            self._model.save_model(str(model_path))
            serialize(self._labels, label_path)
            serialize(self._label_encoder, label_encoder_path)
            params_path.write_text(json.dumps(self._params, indent=4))
        else:
            raise ValueError("Model is not trained or loaded yet.")

    def stop(self) -> None:
        self._gen.stop()

    @staticmethod
    def _serialize_dmatrix(dmatrix: xgb.DMatrix, file_path: str | Path) -> None:
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        dmatrix.save_binary(file_path)

    @staticmethod
    def _deserialize_dmatrix(file_path: str | Path) -> xgb.DMatrix:
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"{file_path} DMatrix not found.")
        return xgb.DMatrix(file_path)

    def _store_batches(
        self, batch_count: int, storage_name: str, return_dmats: bool = False
    ) -> list[xgb.DMatrix] | None:
        storage_dir = self._save_path / "dmats" / storage_name
        dmats = []
        for i in range(batch_count):
            dmat = self._get_next_batch()
            dmat_path = storage_dir / f"{i+1:04}.dmat"

            counter = 1
            while dmat_path.exists():
                # change name if needed
                dmat_path = storage_dir / f"{i+1:04}_{counter:04}.dmat"
                counter += 1

            self._serialize_dmatrix(dmat, dmat_path)
            if return_dmats:
                dmats.append(dmat)
        if return_dmats:
            return dmats

    def _stored_batches_generator(
        self, storage_name: str
    ) -> Generator[xgb.DMatrix, None, None]:
        storage_dir = self._save_path / "dmats" / storage_name
        if not storage_dir.is_dir():
            raise NotADirectoryError(f"{storage_dir} is not a directory.")

        for file_path in storage_dir.iterdir():
            if file_path.is_file():
                dmatrix = self._deserialize_dmatrix(file_path)
                yield dmatrix

    def _get_next_batch(self, storage_name: str | None = None) -> xgb.DMatrix:
        if self._last_batch_id is None:
            self._last_batch_id = -1
        self._last_batch_id += 1
        try:
            if storage_name:
                # store before returning
                return self._store_batches(
                    batch_count=1, storage_name=storage_name, return_dmats=True
                )[0]
            return next(self._gen.get())
        except StopIteration:
            LOG.exception(
                f"No more batch available: {self._last_batch_id + 1} / {self._max_batch_count}"
            )
            raise

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
