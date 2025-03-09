import json
from pathlib import Path
import pickle
from xgboost import XGBClassifier, DMatrix, Booster

from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
import numpy as np
from tqdm.auto import tqdm
from api import API
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
        # self.trials = Trials()


    # def get_model()

    def train(
        self, dtrain: DMatrix, kfold: int | None = 5, tax_id_2_name: bool = True
    ) -> dict | list:
        "When using k-fold, just return metrics,"
        "otherwise, return nothing (and keep model in self._model + labels in self._label)."
        params = self._params.copy()

        X = dtrain.get_data()
        y = dtrain.get_label().astype(int)
        label_encoder = LabelEncoder()
        y_encoded = label_encoder.fit_transform(y)
        labels = label_encoder.classes_
        if tax_id_2_name and self._api:
            labels = [
                f'{self._api[tax_id]["ScientificName"]} [{tax_id}]' for tax_id in labels
            ]

        params["seed"] = 2025
        if self._use_gpu:
            params["tree_method"] = "hist"
            params["device"] = "cuda"
            # can't test it on my laptop(can't install cupy)
            X_dense = X.toarray() if hasattr(X, "toarray") else X
            X = torch.tensor(X_dense, device="cuda")
            y_dense = y.toarray() if hasattr(y, "toarray") else X
            y = torch.tensor(y_dense, device="cuda")

        if kfold is not None:
            kf = KFold(n_splits=kfold, shuffle=True, random_state=2025)
            y_pred = np.zeros(y_encoded.shape)

            for train_index, valid_index in tqdm(kf.split(X), desc="k-fold"):
                model = XGBClassifier(**params)
                model.fit(X[train_index], y_encoded[train_index])
                preds = model.predict(X[valid_index])
                y_pred[valid_index] = preds

            report = classification_report(
                y_encoded, y_pred, target_names=labels, output_dict=True
            )
            return report
        else:
            final_model = XGBClassifier(**params)
            final_model.fit(X, y_encoded)
            self._model = final_model
            self._labels = labels

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
