"""Creates the XGB models"""
import os
from copy import copy
from xgboost import DMatrix, train
from xgboost.core import XGBoostError
import uuid

from wisp_light.dataset.refSeqDataset import TAXO_LEVELS


def make_model(output_dir: str, filtered_reads, mappings_data,
               # database: dict
               params: dict, logger,  taxo_level: str, taxo_target: str) :
    """
    Builds and trains an XGBoost model to classify reads at the next taxonomic level.

    Parameters
    ----------
    output_dir : str
        Directory where the model and configuration files will be saved.
    filtered_reads : iterable
        Iterable of (data_by_genome, read) tuples, where `read` is a dict of k-mer counts.
    mappings_data : dict
        Dictionary mapping taxonomic labels to numerical class indices.
    params : dict
        Dictionary of hyperparameters for the XGBoost model, including 'num_rounds_boosting'.
    logger : Logger
        Logger object used to output warnings and errors.
    taxo_level : str
        Current taxonomic level used to filter data (e.g., 'phylum').
    taxo_target : str
        The specific taxonomic label to train the model on at the current level.

    Returns
    -------
    str or None
        Path to the saved XGBoost model file if training succeeds, otherwise `None`.
    """
    # Creating the booster
    # levels_old = ["root", "domain", "phylum", "group", "order", "family", "specie"]
    levels  = ['root'] + TAXO_LEVELS #, 'phylum', 'class', 'order', 'family']
    level_up: int = levels.index(taxo_level) + 1
    next_level: str = levels[level_up]

    if not taxo_level in mappings_data and not taxo_level == 'root':
        raise ValueError(f"Database does not contain {taxo_level} level.")

    mappings: dict = copy(mappings_data[next_level])
    try:
        number_taxa = mappings.pop('number_taxa')
    except KeyError:
        logger.error("Key 'number_taxa was missing, for unknown reason")
        number_taxa = len(mappings)         # Defining default value


    temp_dir = f"{output_dir}/tmp"
    model_dir = f"{output_dir}/model"
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    # We will be creating temporary LibSVM files in order to make our model learn on those,
    # then destroy files in order to save space

    temp_dataset = f"{temp_dir}/{taxo_target}_{taxo_level}_{uuid.uuid4().hex}.txt"

    with open(temp_dataset, 'w', encoding='utf-8') as libsvm_writer:
        for data_by_genome, read in filtered_reads:
            try:
                label_next_level = data_by_genome[next_level]
                id_label_next_level = mappings[label_next_level]
                kmer_pairs = ' '.join(f"{k}:{v}" for k, v in read.items())
                line = f"{id_label_next_level} {kmer_pairs} #{label_next_level}\n"
                libsvm_writer.write(line)
            except KeyError as e:
                logger.warning(f"Missing mapping for label {e} in {next_level} — skipping read.")

        # for data_by_genome in database['datas']:
        #     if taxo_level == 'root' or data_by_genome[taxo_level] == taxo_target: # aps de root dans data_by_genome
        #
        #         for read in data_by_genome['datas']:
        #             id_label_next_level = mappings[data_by_genome[next_level]]
        #             label_next_level = data_by_genome[next_level]
        #             kmer_pairs = ' '.join([str(k) + ':' + str(v) for k, v in read.items()]) # Each read is a dict with code:count for kmer
        #             line = f"{id_label_next_level} {kmer_pairs} #{label_next_level}\n" # sert pour l'eval et non xgboost
        #             libsvm_writer.write(line)

    model_output_path = f"{model_dir}/{taxo_target}_{taxo_level}.json"
    config_output_path = model_output_path.replace('.json', '_params.json')

    model_params = {key: params[key] for key in
                    {"eval_metric", "tree_method", "device", "booster", "objective",  "eta", "max_depth", "random_state"}
                    if key in params}
    model_params['num_class'] = number_taxa
    model_params['nthread'] = 1

    try:
        # Creating the model
        dm_train = DMatrix(temp_dataset+"?format=libsvm")
        bst = train(model_params,dm_train , params['num_rounds_boosting']) # Booster
        bst.save_model(model_output_path)  # Saving the model and its params        # Must go to model_dir

    except XGBoostError as e:
        logger.error(model_params)
        logger.error(f"Error XgBoost  pour {taxo_target} {e}")    # Invalid dataset, we don't want to keep current level
        return None

    with open(config_output_path, 'w', encoding='utf-8') as jwriter:
        jwriter.write(bst.save_config())

    os.remove(temp_dataset)     # Destroying the temporary directory and its contents
    return model_output_path  # Returning the target file


