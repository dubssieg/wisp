"""Creates the XGB models"""
import os
import logging
from utils import setup_logger
from copy import copy
from xgboost import config_context, DMatrix, train
from xgboost.core import XGBoostError
import uuid

logger = setup_logger(__name__, level=logging.INFO)


def make_model(output_dir: str, datas: dict, params: dict, taxo_level: str, taxo_target: str,) :
    """Builds the model and saves it"""
    # Creating the booster
    config_context(booster='gbtree', min_child_weight=1, tree_method='approx', predictor='cpu_predictor')
    levels = ["root", "domain", "phylum", "group", "order", "family", "specie"]
    level_up: int = levels.index(taxo_level) + 1
    next_level: str = levels[level_up]

    if not taxo_level in datas['mappings'] and not taxo_level == 'root':
        raise ValueError(f"Database does not contain {taxo_level} level.")

    mappings: dict = copy(datas['mappings'][next_level])
    try:
        number_taxa: int = mappings.pop('number_taxa')
    except KeyError:
        logger.error("Key 'number_taxa was missing, for unknown reason")
        number_taxa: int = len(mappings)         # Defining default value


    temp_dir = f"{output_dir}/tmp"
    model_dir = f"{output_dir}/model"
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    # We will be creating temporary LibSVM files in order to make our model learn on those,
    # then destroy files in order to save space

    temp_dataset = f"{temp_dir}/{taxo_target}_{taxo_level}_{uuid.uuid4().hex}.txt"

    with open(temp_dataset, 'w', encoding='utf-8') as libsvm_writer:
        for data_by_genome in datas['datas']:
            if taxo_level == 'root' or data_by_genome[taxo_level] == taxo_target: # aps de root dans data_by_genome
                # Sample should be kept for model
                for read in data_by_genome['datas']:
                    id_label_next_level = mappings[data_by_genome[next_level]]
                    label_next_level = data_by_genome[next_level]
                    kmer_pairs = ' '.join([str(k) + ':' + str(v) for k, v in read.items()]) # Each read is a dict with code:count for kmer
                    line = f"{id_label_next_level} {kmer_pairs} #{label_next_level}\n" # sert pour l'eval et non xgboost
                    libsvm_writer.write(line)

    model_output_path = f"{model_dir}/{taxo_target}_{taxo_level}.json"
    config_output_path = model_output_path.replace('.json', '_params.json')

    model_params = {key: params[key] for key in
                    {"eval_metric", "tree_method", "device", "booster", "objective",  "eta", "max_depth"}
                    if key in params}
    model_params['num_class'] = number_taxa

    with open(temp_dataset, 'r', encoding='utf-8') as f:
        content = f.readlines()
    if len(content) == 0:
        logger.error(f"❌ Dataset vide pour {taxo_target} ({taxo_level})")

    try:
        # Creating the model
        bst = train(model_params, DMatrix(temp_dataset+"?format=libsvm"), params['num_rounds_boosting']) # Booster
        bst.save_model(model_output_path)  # Saving the model and its params        # Must go to model_dir

    except XGBoostError as e:
        logger.info(model_params)
        logger.error(f"Error XgBoost  pour {taxo_target} {e}")    # Invalid dataset, we don't want to keep current level
        return None

    with open(config_output_path, 'w', encoding='utf-8') as jwriter:
        jwriter.write(bst.save_config())

    os.remove(temp_dataset)     # Destroying the temporary directory and its contents
    return model_output_path  # Returning the target file


