from xgboost import XGBClassifier, DMatrix
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
import numpy as np


class XGBoostModel:
    def __init__(self, params=None):
        self._params = params if params is not None else {}
        self.model = None
        self.trials = Trials()

    def train_with_kfold(
        self,
        dtrain: DMatrix,
        n_estimators: int = 100,
        nfold: int = 5,
    ) -> dict:
        X = dtrain.get_data()
        y = dtrain.get_label().astype(int)

        label_encoder = LabelEncoder()
        y_encoded = label_encoder.fit_transform(y)
        labels = label_encoder.classes_

        kf = KFold(n_splits=nfold, shuffle=True, random_state=2025)
        y_pred = np.zeros(y_encoded.shape)

        for train_index, valid_index in kf.split(X):
            model = XGBClassifier(**self._params, n_estimators=n_estimators)
            model.fit(X[train_index], y_encoded[train_index])
            preds = model.predict(X[valid_index])
            y_pred[valid_index] = preds

        report = classification_report(
            y_encoded, y_pred, target_names=labels, output_dict=True
        )

        final_model = XGBClassifier(**self._params, n_estimators=n_estimators)
        final_model.fit(X, y_encoded)

        return report
