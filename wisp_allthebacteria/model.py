from datetime import datetime
import json
from pathlib import Path
import pickle
import time
from typing import Generator
from matplotlib import pyplot as plt
import seaborn as sns
from database import DMatrixGeneratorFactory, Database

# from xgboost import DMatrix, Booster
import xgboost as xgb

from sklearn.model_selection import KFold

# from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix

# from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
import numpy as np
from tqdm.auto import tqdm
from api import API
from utils import system_stats, format_duration, format_size

# import torch

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


class XGBoostModel:
    def __init__(
        self,
        params=None,
        api: API | None = None,
        use_gpu: bool = True,
        scientific_name: bool = True,
    ):
        self._params = params if params is not None else DEFAULT_PARAMETERS
        self._api = api
        self._use_gpu = use_gpu
        self._model = None
        self._tax_id_labels = None
        self._labels = None
        self._report = None
        self._scientific_name = scientific_name

        # self.trials = Trials()

    def train(
        self,
        # dtrain: xgb.DMatrix | Generator[xgb.DMatrix, None, None],
        dtrain_factory: DMatrixGeneratorFactory,
        tax_id_classes: list[int],
        num_boost_round: int,
        kfold: int | None = None,
    ) -> dict:
        "Train (full data) or evaluate (kfold) and return report. Keep model, labels, report."
        params = self._params.copy()
        if "seed" not in params:
            params["seed"] = 2025
        params["num_class"] = len(tax_id_classes)

        self._report = dict()
        self._report["start_dt"] = datetime.now()
        total_duration_start = time.time()

        label_encoder = self._update_labels(tax_id_classes)

        self._report["params"] = params
        self._report["num_rows"] = 0
        self._report["x_size"] = 0
        self._report["num_boost_round"] = num_boost_round
        self._report["labels"] = self._labels

        if kfold is None:
            self._full_train(
                dtrain=dtrain_factory.make_dmatrix_generator(),
                num_boost_round=num_boost_round,
                label_encoder=label_encoder,
                params=params,
            )
        else:
            pass

        self._report["end_dt"] = datetime.now()
        self._report["total_duration"] = time.time() - total_duration_start

        # for each batch

        # can't test since I can't install cupy TODO: try something else
        # if self._use_gpu:
        #     params["tree_method"] = "hist"
        #     params["device"] = "cuda"

        #     X_dense = X.toarray() if hasattr(X, "toarray") else X
        #     X = torch.tensor(X_dense, device="cuda")
        #     y_dense = y.toarray() if hasattr(y, "toarray") else X
        #     y = torch.tensor(y_dense, device="cuda")

        # need major refactoring for generator

        # fold_start_time = time.time()
        # kf = KFold(n_splits=kfold, shuffle=True, random_state=2025)
        # y_pred = np.zeros(y_encoded.shape)

        # for train_index, valid_index in tqdm(
        #     kf.split(X), desc="k-fold", total=kfold
        # ):
        #     model = XGBClassifier(**params)
        #     model.fit(X[train_index], y_encoded[train_index])
        #     cpu_mem_stats.append(system_stats())
        #     y_pred[valid_index] = model.predict(X[valid_index])
        #     train_durations.append(time.time() - fold_start_time)

        # self._report["classification_report"] = classification_report(
        #     y_encoded, y_pred, target_names=self._labels, output_dict=True
        # )
        # self._report["confusion_matrix"] = confusion_matrix(
        #     y_encoded, y_pred, labels=range(len(self._labels))
        # )

        return self._report

    def _full_train(
        self,
        dtrain: xgb.DMatrix,
        params: dict,
        num_boost_round: int,
        label_encoder: dict,
    ) -> dict:
        cpu_mem_stats = list()
        # if dtrain in not a generator, make it look like one
        if isinstance(dtrain, xgb.DMatrix):
            dtrain = [dtrain]

        # XGBClassifier do not fully support incremental training
        # self._model = XGBClassifier(**params)
        self._model = None
        for i, dtrain_i in enumerate(dtrain):

            if i == 0:
                self._report["train_dtype"] = dtrain_i.get_data().dtype
            self._report["num_rows"] += dtrain_i.num_row()
            self._report["x_size"] += (
                dtrain_i.num_row()
                * dtrain_i.num_col()
                * dtrain_i.get_data().dtype.itemsize
            )
            y = dtrain_i.get_label().astype(int)
            y_encoded = [label_encoder[label] for label in y]
            dtrain_i_encoded = xgb.DMatrix(dtrain_i.get_data(), label=y_encoded)

            self._model = xgb.train(
                params,
                dtrain_i_encoded,
                num_boost_round=num_boost_round,
                xgb_model=self._model,
            )
            cpu_mem_stats.append(system_stats())

        self._report["system_stats"] = cpu_mem_stats

    def _update_labels(self, tax_id_classes: list[int]) -> dict:
        self._tax_id_labels = tax_id_classes
        if self._scientific_name and self._api:
            self._labels = [
                f'{self._api[tax_id]["ScientificName"]} [{tax_id}]'
                for tax_id in self._tax_id_labels
            ]

        return {label: idx for idx, label in enumerate(self._tax_id_labels)}

    # def _update_labels(self, y: np.array) -> dict:
    #     unique_tax_id_labels = np.unique(y)

    #     for tax_id_label in unique_tax_id_labels:
    #         if tax_id_label not in self._tax_id_labels:
    #             self._tax_id_labels.append(tax_id_label)

    #     if self._get_scientific_name and self._api:
    #         self._labels = [
    #             f'{self._api[tax_id]["ScientificName"]} [{tax_id}]'
    #             for tax_id in self._tax_id_labels
    #         ]
    #     else:
    #         self._labels = self._tax_id_labels

    #     return {label: idx for idx, label in enumerate(self._tax_id_labels)}

    def save_report(self, additional_data: dict, dir_path: str | Path):
        """save self._report + additional data"""
        report = {**self._report, **additional_data}

        dt_format = "%Y-%m-%d %H:%M:%S"
        report_lines = []

        # header
        if "header" in report:
            for k, v in report["header"].items():
                report_lines.append(f"{k} : {v}")
            report_lines.append("")

        # basic
        report_lines.append("=== Rapport d'entraînement / évaluation ===")
        report_lines.append(f"Date de début : {report['start_dt'].strftime(dt_format)}")
        report_lines.append(f"Date de fin : {report['end_dt'].strftime(dt_format)}")
        report_lines.append(
            f"Durée totale : {format_duration(report['total_duration'])}"
        )
        report_lines.append("")

        report_lines.append("\n=== Paramètres ===")
        report_lines.append(f"Itérations de boosting : {report['num_boost_round']}")
        report_lines.append(json.dumps(report["params"], indent=4))
        report_lines.append("")

        # training
        report_lines.append("\n=== Données d'entraînement ===")
        report_lines.append(f"Nombre de séquences : {report['num_rows']}")
        report_lines.append(f"Type de données : {report['train_dtype']}")
        report_lines.append(f"Taille estimée : {format_size(report['x_size'])}")
        report_lines.append("")

        # classification
        report_lines.append("\n=== Rapport de classification ===")
        if "classification_report" in report:

            for label, metrics in report["classification_report"].items():
                if isinstance(metrics, dict):
                    report_lines.append(f"Label : {label}")
                    for metric_name, value in metrics.items():
                        report_lines.append(f"  {metric_name} : {value:.4f}")
                else:
                    report_lines.append(f"{label} : {metrics:.4f}")
        else:
            report_lines.append("Pas de rapport de classification disponible.")
        report_lines.append("")

        # confusion matrix
        report_lines.append("\n=== Matrice de confusion ===")
        if "confusion_matrix" in report:
            report_lines.append(
                self._confusion_matrix_to_ascii(
                    report["confusion_matrix"], report["labels"]
                )
            )
        else:
            report_lines.append("Pas de matrice de confusion disponible.")
        report_lines.append("")

        # system
        report_lines.append("\n=== Statistiques Système ===")
        for i, stats in enumerate(report["system_stats"]):
            if len(report["system_stats"]) > 1:
                report_lines.append(f" - Batch {i + 1}:")
            report_lines.append(f"  Utilisation CPU : {stats['cpu_usage']}%")
            report_lines.append(f"  Utilisation RAM : {stats['ram_usage_percent']}%")
            report_lines.append(f"  RAM totale : {format_size(stats['total_ram'])}")
            report_lines.append(f"  RAM utilisée : {format_size(stats['used_ram'])}")
            report_lines.append("")

        # write
        dir_path = Path(dir_path)
        dir_path.mkdir(parents=True, exist_ok=True)
        txt_path = dir_path / "report.txt"
        cm_path = dir_path / "confusion_matrix.png"

        with open(txt_path, "w") as file:
            file.write("\n".join(report_lines))

        if "confusion_matrix" in report:
            self._plot_and_save_confusion_matrix(
                report["confusion_matrix"], report["labels"], cm_path
            )

    @staticmethod
    def _confusion_matrix_to_ascii(conf_matrix: list, labels: list) -> str:
        max_label_length = max(len(label) for label in labels)
        header = f"{'':>{max_label_length+2}} " + " ".join(
            f"{label:>{max_label_length}}" for label in labels
        )
        lines = [header]

        for i, label in enumerate(labels):
            line = f"{label:{max_label_length}}: " + " ".join(
                f"{conf_matrix[i, j]:{max_label_length}}" for j in range(len(labels))
            )
            lines.append(line)
        return "\n".join(lines)

    @staticmethod
    def _plot_and_save_confusion_matrix(
        conf_matrix: list, labels: list, cm_path: str | Path
    ) -> None:
        _, ax = plt.subplots(figsize=(10, 7))
        sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", ax=ax, cbar=False)

        ax.set_xlabel("Predicted Labels")
        ax.set_ylabel("True Labels")
        ax.set_title("Confusion Matrix")

        ax.set_xticks(np.arange(len(labels)) + 0.5)
        ax.set_yticks(np.arange(len(labels)) + 0.5)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels, rotation=0)

        plt.tight_layout()
        plt.savefig(cm_path)
        plt.close()

    def load(self, dir_path: str | Path) -> None:
        """Load model from file."""
        model_path = Path(dir_path) / "model.bin"
        label_path = Path(dir_path) / "labels.pkl"
        params_path = dir_path / "params.json"
        if not model_path.is_file():
            raise FileNotFoundError(f"Model file not found: {model_path}")
        if not label_path.is_file():
            raise FileNotFoundError(f"Label file not found: {label_path}")
        if not params_path.is_file():
            raise FileNotFoundError(f"Param file not found: {label_path}")

        self._model = xgb.Booster()
        self._model.load_model(str(model_path))
        with open(label_path, "rb") as f:
            self._model = pickle.load(f)
        self._params = json.loads(params_path.read_text())

    def save(self, dir_path: str | Path) -> None:
        """Save model to directory."""
        dir_path = Path(dir_path)
        dir_path.mkdir(parents=True, exist_ok=True)
        model_path = dir_path / "model.bin"
        label_path = dir_path / "labels.pkl"
        params_path = dir_path / "params.json"

        if self._model is not None:
            self._model.save_model(str(model_path))
            with open(label_path, "wb") as f:
                pickle.dump(self._labels, f)
            params_path.write_text(json.dumps(self._params, indent=4))
        else:
            raise ValueError("Model is not trained or loaded yet.")

    def optimize_hyperparameters(self, dtrain: xgb.DMatrix, nfold=5):
        """Optimize hyperparameters using Hyperopt."""
        # TODO do it again

    def predict(self, dtest: xgb.DMatrix):
        """Make predictions on the test set."""
        if self._model is None:
            raise ValueError("Model loaded.")
        return self._model.predict(dtest)
