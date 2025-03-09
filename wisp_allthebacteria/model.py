from datetime import datetime
import json
from pathlib import Path
import pickle
import time
from matplotlib import pyplot as plt
import seaborn as sns
from xgboost import XGBClassifier, DMatrix, Booster

from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
import numpy as np
from tqdm.auto import tqdm
from api import API
from utils import system_stats, format_duration, format_size
import torch

DEFAULT_PARAMETERS = {
    "n_estimators": 100,  # 1000 + early_stopping ?
    "max_bin": 256,  # for gpu_hist
    "grow_policy": "depthwise",  # better with big dataset
    "objective": "multi:softmax",
    "eval_metric": "logloss",
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
    def __init__(self, params=None, api: API | None = None, use_gpu: bool = True):
        self._params = params if params is not None else DEFAULT_PARAMETERS
        self._api = api
        self._use_gpu = use_gpu
        self._model = None
        self._labels = None
        self._report = None

        # self.trials = Trials()

    # def get_model()

    def train(
        self, dtrain: DMatrix, kfold: int | None = None, tax_id_2_name: bool = True
    ) -> dict:
        "Train (full data) or evaluate (kfold) and return report. Keep model, labels, report."
        params = self._params.copy()
        self._report = dict()
        self._report["start_dt"] = datetime.now()
        total_duration_start = time.time()

        X = dtrain.get_data()
        y = dtrain.get_label().astype(int)
        label_encoder = LabelEncoder()
        y_encoded = label_encoder.fit_transform(y)
        self._labels = label_encoder.classes_
        if tax_id_2_name and self._api:
            self._labels = [
                f'{self._api[tax_id]["ScientificName"]} [{tax_id}]'
                for tax_id in self._labels
            ]

        if "seed" not in params:
            params["seed"] = 2025
        if self._use_gpu:
            params["tree_method"] = "hist"
            params["device"] = "cuda"
            # can't test it on my laptop(can't install cupy)
            X_dense = X.toarray() if hasattr(X, "toarray") else X
            X = torch.tensor(X_dense, device="cuda")
            y_dense = y.toarray() if hasattr(y, "toarray") else X
            y = torch.tensor(y_dense, device="cuda")

        train_durations = []
        cpu_mem_stats = []

        self._report["params"] = params
        self._report["train_shape"] = X.shape
        self._report["train_dtype"] = X.dtype
        self._report["train_estimated_size_octets"] = (
            X.shape[0] * X.shape[1] * X.dtype.itemsize
        )

        if kfold is not None:
            fold_start_time = time.time()
            kf = KFold(n_splits=kfold, shuffle=True, random_state=2025)
            y_pred = np.zeros(y_encoded.shape)

            for train_index, valid_index in tqdm(
                kf.split(X), desc="k-fold", total=kfold
            ):
                model = XGBClassifier(**params)
                model.fit(X[train_index], y_encoded[train_index])
                cpu_mem_stats.append(system_stats())
                y_pred[valid_index] = model.predict(X[valid_index])
                train_durations.append(time.time() - fold_start_time)

            self._report["classification_report"] = classification_report(
                y_encoded, y_pred, target_names=self._labels, output_dict=True
            )
            self._report["confusion_matrix"] = confusion_matrix(
                y_encoded, y_pred, labels=range(len(self._labels))
            )

        else:
            train_start_time = time.time()
            final_model = XGBClassifier(**params)
            final_model.fit(X, y_encoded)
            cpu_mem_stats.append(system_stats())
            self._model = final_model
            train_durations.append(time.time() - train_start_time)

        self._report["labels"] = self._labels
        self._report["system_stats"] = cpu_mem_stats
        self._report["end_dt"] = datetime.now()
        self._report["total_duration"] = time.time() - total_duration_start

        return self._report

    def save_report(self, additional_data: dict, dir_path: str | Path):
        """save self._report + additional data"""
        report = {**self._report, **additional_data}

        dt_format = "%Y-%m-%d %H:%M:%S"
        report_lines = []

        # header
        if "header" in report:
            for k, v in report["header"].items():
                report_lines.append(f"{k} : {v}")

        # basic
        report_lines.append("=== Rapport d'entraînement / évaluation ===")
        report_lines.append(f"Date de début : {report['start_dt'].strftime(dt_format)}")
        report_lines.append(f"Date de fin : {report['end_dt'].strftime(dt_format)}")
        report_lines.append(
            f"Durée totale : {format_duration(report['total_duration'])}"
        )
        report_lines.append("")

        report_lines.append("\n=== Paramètres ===")
        report_lines.append(json.dumps(report["params"], indent=4))
        report_lines.append("")

        # training
        report_lines.append("\n=== Données d'entraînement ===")
        report_lines.append(f"Format des données : {report['train_shape']}")
        report_lines.append(f"Type de données : {report['train_dtype']}")
        report_lines.append(
            f"Taille estimée : {format_size(report['train_estimated_size_octets'])}"
        )
        report_lines.append("")

        # classification
        report_lines.append("\n=== Rapport de classification ===")
        for label, metrics in report["classification_report"].items():
            if isinstance(metrics, dict):
                report_lines.append(f"Label : {label}")
                for metric_name, value in metrics.items():
                    report_lines.append(f"  {metric_name} : {value:.4f}")
            else:
                report_lines.append(f"{label} : {metrics:.4f}")
        report_lines.append("")

        # confusion matrix
        report_lines.append("\n=== Matrice de confusion ===")
        report_lines.append(
            self._confusion_matrix_to_ascii(
                report["confusion_matrix"], report["labels"]
            )
        )
        report_lines.append("")

        # system
        report_lines.append("\n=== Statistiques Système ===")
        for index, stats in enumerate(report["system_stats"]):
            if len(report["system_stats"]) > 1:
                report_lines.append(f" - Fold {index + 1}:")
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

        self._model = Booster()
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

    def optimize_hyperparameters(self, dtrain: DMatrix, nfold=5):
        """Optimize hyperparameters using Hyperopt."""
        # TODO do it again

    def predict(self, dtest: DMatrix):
        """Make predictions on the test set."""
        if self._model is None:
            raise ValueError("Model loaded.")
        return self._model.predict(dtest)
