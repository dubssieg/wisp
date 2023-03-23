"Creates the XGB models"
from os import path
from pathlib import Path
from typing import Any
from xgboost import Booster, config_context, DMatrix, train
# from shutil import rmtree


def make_model(
        datas: Any,
        model_name: str,
        classification_level: str,
        target_dataset: str,
        num_rounds_boosting: int = 10,
        eta: float = 0.3,
        maximum_depth: int = 10
) -> tuple[str, str]:
    "Builds the model and saves it"
    # Creating the booster
    config_context(
        booster='gbtree',
        min_child_weight=1,
        tree_method='approx',
        predictor='cpu_predictor'
    )

    next_level: str = (levels := ["root", "domain", "phylum", "group", "order", "family", "specie"])[
        (levels).index(classification_level)+1]

    if not classification_level in datas['mappings'] and classification_level != 'root':
        raise ValueError(
            f"Database does not contain {classification_level} level.")

    mappings: dict = datas['mappings'][next_level]
    number_taxa: int = mappings.pop('number_taxa')

    model_parameters: dict = {
        'max_depth': maximum_depth,
        'objective': 'multi:softprob',
        'num_class': number_taxa,
        'eta': eta,
        'eval_metric': 'mlogloss'
    }

    # We will be creating temporary LibSVM files in order to make our model learn on those,
    # then destroy files in order to save space
    Path(temp_dir :=
         f"{path.dirname(__file__)}/tmp").mkdir(parents=True, exist_ok=True)

    # We create the dir to store the database
    Path(model_dir := f"{path.dirname(__file__)}/model/{Path(model_name).stem}").mkdir(
        parents=True, exist_ok=True)

    # We create temporary files
    with open(temp_dataset := f"{temp_dir}/{target_dataset}.txt", 'w', encoding='utf-8') as libsvm_writer:
        for sample in datas['datas']:
            if classification_level == 'root' or sample[classification_level] == target_dataset:
                # Sample should be kept for model
                for read in sample['datas']:
                    # Each read is a dict with code:count for kmer
                    libsvm_writer.write(
                        f"{mappings[sample[next_level]]} {' '.join([str(k)+':'+str(v) for k,v in read.items()])} #{sample[next_level]}")

    # Creating the model
    bst: Booster = train(
        model_parameters,
        DMatrix(temp_dataset),
        num_rounds_boosting
    )

    # Saving the model and its params
    # Must go to model_dir
    bst.save_model(model_output_path := f"{model_dir}/{target_dataset}.json")

    with open(config_output_path := f"{model_dir}/{target_dataset}_params.json", 'w', encoding='utf-8') as jwriter:
        jwriter.write(bst.save_config())

    # Destroying the temporary directory and its contents
    # rmtree(temp_dir)

    # returning the target file
    return model_output_path, config_output_path
