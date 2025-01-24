"""Creates the XGB models"""
import os
from typing import Any
from copy import copy
from xgboost import Booster, config_context, DMatrix, train
from xgboost.core import XGBoostError


def make_model(output_dir: str,
                datas: Any,
                model_name: str,
                classification_level: str,
                target_dataset: str,
                num_rounds_boosting: int = 10,
                eta: float = 0.3,
                maximum_depth: int = 10,) :
    """Builds the model and saves it"""
    # Creating the booster
    config_context(booster='gbtree', min_child_weight=1, tree_method='approx', predictor='cpu_predictor')
    levels = ["root", "domain", "phylum", "group", "order", "family", "specie"]
    level_up = levels.index(classification_level)+1
    next_level: str = levels[level_up]

    if not classification_level in datas['mappings'] and not classification_level == 'root':
        raise ValueError(f"Database does not contain {classification_level} level.")

    mappings: dict = copy(datas['mappings'][next_level])
    try:
        number_taxa: int = mappings.pop('number_taxa')
    except KeyError:
        print("Key 'number_taxa was missing, for unknown reason")
        number_taxa: int = len(mappings)         # Defining default value


    model_parameters: dict = {'max_depth': maximum_depth, 'objective': 'multi:softprob', 'num_class': number_taxa,
                              'eta': eta,'eval_metric': 'mlogloss' }
    try:
        # We will be creating temporary LibSVM files in order to make our model learn on those,
        # then destroy files in order to save space
        temp_dir = f"{output_dir}/tmp"
        model_dir = f"{output_dir}/model/{os.path.splitext(os.path.basename(model_name))[0]}"
        os.makedirs(model_dir, exist_ok=True)
        os.makedirs(temp_dir, exist_ok=True)

        temp_dataset = f"{temp_dir}/{target_dataset}_{classification_level}.txt"

        with open(temp_dataset, 'w', encoding='utf-8') as libsvm_writer:
            for sample in datas['datas']:
                if classification_level == 'root' or sample[classification_level] == target_dataset:
                    # Sample should be kept for model
                    for read in sample['datas']:
                        # Each read is a dict with code:count for kmer
                        libsvm_writer.write(f"{mappings[sample[next_level]]} {' '.join([str(k)+':'+str(v) for k,v in read.items()])} #{sample[next_level]}\n")
        # Creating the model
        bst = train(model_parameters, DMatrix(temp_dataset+"?format=libsvm"), num_rounds_boosting) # Booster

        # Saving the model and its params        # Must go to model_dir
        model_output_path = f"{model_dir}/{target_dataset}_{classification_level}.json"
        bst.save_model(model_output_path)
        config_output_path = f"{model_dir}/{target_dataset}_{classification_level}_params.json"

        with open(config_output_path, 'w', encoding='utf-8') as jwriter:
            jwriter.write(bst.save_config())

        os.remove(temp_dataset)     # Destroying the temporary directory and its contents
        return model_output_path, config_output_path  # Returning the target file

    except XGBoostError as e:
        # Invalid dataset, we don't want to keep current level
        print("Error XgBoost ", e)
        return None, None
