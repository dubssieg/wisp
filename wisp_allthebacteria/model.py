import json
import logging
from pathlib import Path
from typing import Generator
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
        save_path: str | Path,
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
        self._save_path = save_path = Path(save_path).resolve()
        self._model = None
        self._labels = None
        self._gen = None
        self._last_batch_id = None
        self._max_batch_count = None
        self._label_encoder = LabelEncoder()
        LOG.debug(
            f"XGBoostModel({rank=}, {batch_size=}, {normalize=}, {seed=}, {generator_threads=}, {save_path=})"
        )

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
        self._max_batch_count = (
            self._gen.available_batches_count()
        )  # TODO: rework after balancing options
        if batch_count is None:
            batch_count = self._max_batch_count

        if batch_count > self._max_batch_count:
            raise ValueError(
                f"not enough batches available: {batch_count} > {self._max_batch_count}"
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
            dtrain = self._get_next_batch()
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

    # def evaluate_dmatrix(self, dtest: xgb.DMatrix) -> dict:

    #     y_true = dtest.get_label().astype(int)
    #     y_true_encoded = self._label_encoder.transform(y_true)
    #     dtest_encoded = xgb.DMatrix(dtest.get_data(), label=y_true_encoded)
    #     y_pred_encoded = self._model.predict(dtest_encoded).astype(int)

    def evaluate(self, batch_count: int, use_stored_batches: str | None = None) -> dict:
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
            dtest = self._get_next_batch(storage_name="test")
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

        LOG.debug("Model evaluated")
        return self._report

    def load(self) -> None:
        """Load model from file."""
        model_path = self._save_path / "model.bin"
        label_path = self._save_path / "labels.pkl"
        label_encoder_path = self._save_path / "labels_encoder.pkl"
        params_path = self._save_path / "params.json"
        for f in (model_path, label_path, label_encoder_path, params_path):
            if not f.is_file:
                raise FileNotFoundError(f)

        self._model = xgb.Booster()
        self._model.load_model(str(model_path))
        self._labels = deserialize(label_path)
        self._params = json.loads(params_path.read_text())

    def save(self) -> None:
        """Save model to directory."""
        self._save_path.mkdir(parents=True, exist_ok=True)
        model_path = self._save_path / "model.bin"
        label_path = self._save_path / "labels.pkl"
        label_encoder_path = self._save_path / "labels_encoder.pkl"
        params_path = self._save_path / "params.json"

        if self._model is not None:
            self._model.save_model(str(model_path))
            serialize(self._labels, label_path)
            serialize(self._labels_encoder, label_encoder_path)
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
                f"No more batch available: {self._last_batch_id} / {self._max_batch_count}"
            )
            raise
