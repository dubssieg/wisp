import json
import logging
from pathlib import Path
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder
import xgboost as xgb
from utils import serialize, deserialize

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


class XGBoostModel:
    def __init__(
        self,
        rank: str,
        database: Database,
        params: dict | None = None,
        normalize: str | None = None,
        batch_size: int = 1_000_000,
        seed: int = 2025,
        api: API | None = None,
        generator_threads: int = 10,
    ):
        self._params = params if params is not None else DEFAULT_PARAMETERS
        self._seed = seed
        self._rank = rank
        self._normalize = normalize
        self._batch_size = batch_size
        self._api = api
        self._database = database
        self._generator_threads = generator_threads
        self._model = None
        self._labels = None
        self._gen = None
        self._last_batch_id = None
        self._label_encoder = LabelEncoder()

    def train(
        self,
        batch_count: int | None = None,
        num_boost_round: int = 100,
    ) -> None:
        """Train a model"""

        LOG.debug("Training model...")

        # sample generator
        self._gen = ByRankGenerator(
            database=self._database,
            api=self._api,
            rank=self._rank,
            batch_size=self._batch_size,
            normalize=self._normalize,
            seed=self._seed,
            buffer_threads=self._generator_threads,
        )
        max_batch_count = self._gen.available_batches_count()
        if batch_count is None:
            batch_count = max_batch_count

        if batch_count > max_batch_count:
            raise ValueError(
                f"not enough batches available: {batch_count} > {max_batch_count}"
            )

        # labels
        self._labels = self._gen.labels(self._rank)

        self._label_encoder.fit(self._labels)
        # label_encoder.inverse_transform(encoded_labels)

        # params
        params = self._params.copy()
        if "seed" not in params:
            params["seed"] = self._seed
        params["num_class"] = len(self._labels)

        self._gen.start()

        # model
        self._model = None
        for i in range(batch_count):
            self._last_batch_id = i
            dtrain = next(self._gen.get())
            LOG.debug(f"Train batch {i + 1} / {batch_count}")

            y = dtrain.get_label().astype(int)
            y_encoded = self._label_encoder.transform(y)
            dtrain_encoded = xgb.DMatrix(dtrain.get_data(), label=y_encoded)
            self._model = xgb.train(
                params,
                dtrain_encoded,
                num_boost_round=num_boost_round,
                xgb_model=self._model,
            )
        LOG.debug("Model trained")

    def evaluate(self, batch_count: int) -> dict:
        """Evaluate the model on the next available batches."""

        LOG.debug("Evaluating model...")
        if self._model is None:
            raise ValueError("Model has not been trained yet.")

        max_batch_count = self._gen.available_batches_count()
        if self._last_batch_id + batch_count > max_batch_count:
            raise ValueError(
                f"Not enough batches available for evaluation: "
                f"{self._last_batch_id + batch_count} > {max_batch_count}"
            )

        self._report = {}
        all_y_true_encoded = []
        all_y_pred_encoded = []

        # evaluate
        for i in range(batch_count):
            self._last_batch_id += 1
            dtest = next(self._gen.get())
            LOG.debug(
                f"Evaluation batch {i + 1} / {batch_count} (total with training: {self._last_batch_id} / {max_batch_count})"
            )

            y_true = dtest.get_label().astype(int)
            y_true_encoded = self._label_encoder.transform(y_true)
            dtest_encoded = xgb.DMatrix(dtest.get_data(), label=y_true_encoded)
            y_pred_encoded = self._model.predict(dtest_encoded).astype(int)

            all_y_true_encoded.extend(y_true_encoded)
            all_y_pred_encoded.extend(y_pred_encoded)

        y_true = self._label_encoder.inverse_transform(all_y_true_encoded)
        y_pred = self._label_encoder.inverse_transform(all_y_pred_encoded)
        # scientific names
        if self._api:
            sn_map = {
                tax_id: self._api[tax_id]["ScientificName"] for tax_id in self._labels
            }
            y_true = [sn_map[tax_id] for tax_id in y_true]
            y_pred = [sn_map[tax_id] for tax_id in y_pred]

        LOG.debug("Model evaluated, getting report...")

        # generate classification report and confusion matrix
        self._report["classification_report"] = classification_report(
            y_true,
            y_pred,
            output_dict=True,
        )
        self._report["confusion_matrix"] = confusion_matrix(y_true, y_pred)

        LOG.debug("Model report DONE")
        return self._report

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
        self._labels = deserialize(label_path)
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
            serialize(self._labels, label_path)
            params_path.write_text(json.dumps(self._params, indent=4))
        else:
            raise ValueError("Model is not trained or loaded yet.")

    def stop(self) -> None:
        self._gen.stop()
