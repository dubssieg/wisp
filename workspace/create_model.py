from xgboost import Booster, config_context, DMatrix
from json import load


def make_model(
        path_to_database: str,
        classification_level: str,
        target_dataset: str
) -> None:
    "Builds the model and saves it"
    # Creating the booster
    bst: Booster = Booster()
    config_context(
        booster='gbtree',
        min_child_weight=1,
        tree_method='approx',
        predictor='cpu_predictor'
    )
    dtrain: DMatrix = DMatrix()

    # Loading data
    with open(path_to_database, 'r', encoding='utf-8') as jdb:
        datas = load(jdb)

    if not classification_level in datas['mappings']:
        raise ValueError(
            f"Database does not contain {classification_level} level.")

    mappings: dict = datas['mappings'][classification_level]
    number_taxa: int = mappings.pop('number_taxa')

    for sample in datas['datas']:
        if sample[classification_level] == target_dataset:
            # Sample should be kept for model
            pass
