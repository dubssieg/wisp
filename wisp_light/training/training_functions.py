"""Builds predictions from reads"""
import logging
import os
import yaml
import time
import sys
import pickle
import json
from tqdm import tqdm
from treelib import Tree
from Bio import SeqIO

from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed

from create_model import make_model
from create_database import check_parameters
from create_prediction import prediction
from utils import  extract_majority_classification, setup_logger
from metrics import ConfusionMatrixTracker, compute_accuracy_from_conf_matrix_df

sys.path.append('../../..')
from wisp.wisp_light.dataset.bactero_set import TAXO_LEVELS
from wisp.wisp_light.visu.plots_tools import plot_conf_mat


def train_model_targets(phylo_tree, exp_dir, params, logger, num_processes=4):
    start_model = time.time()

    levels = ['root', 'domain', 'phylum', 'group', 'order']

    nodes_per_level: dict = {
        level: [node.tag for node in list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))]
        for i, level in enumerate(levels)}  # jusqu'à order (drop family level)

    # display_json_preview(output_file, num_elements=1)
    with open(f'{exp_dir}/databases.json', 'r', encoding='utf-8') as jdb:
        database = json.load(jdb)  # Loading data => should be put in the main call to escape loading it at each iteration

    logger.info("Starting model creation")
    # datas = {'datas': list_59_data,'mappings': taxa_code_by_level}
    classif_targets = [(taxo_level, taxo_target)
                       for taxo_level, targets in nodes_per_level.items()
                       for taxo_target in targets
                       ]

    make_model_partial = partial(make_model, exp_dir, database, params, logger)
    logger.info(f"Lancement de {len(classif_targets)} modèles avec num_processes={num_processes}")

    with ThreadPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(make_model_partial, *classif_target): classif_target for classif_target in
                   classif_targets}

        for idx, future in enumerate(as_completed(futures)):  # Gestion des tâches dès qu'elles terminent
            taxo_level, taxo_target = futures[future]

            try:
                model_path = future.result()  # Récupération du chemin du modèle généré

                if model_path is not None:
                    key = f"{taxo_target.lower()}_{taxo_level}"
                    try:
                        node = phylo_tree[key]
                        node.data.model_path = os.path.basename(model_path)  # Stockage du chemin relatif
                        logger.info(f"✅ [{idx}/{len(futures)}] Modèle   pour le niveau {taxo_level}: {taxo_target} ")

                    except KeyError:
                        logger.error(f"❌ [{idx}/{len(futures)}] Clé manquante dans l'arbre phylogénétique : {key}."
                                     f" Suppression du nœud {taxo_target.lower()}")
                        phylo_tree.remove_node(taxo_target.lower())

                else:
                    logger.warning(f"⚠ [{idx}/{len(futures)}] Modèle non généré , model_path=None {taxo_target} ({taxo_level})")

            except Exception as e:
                logger.error(f"🔥 Erreur lors de l'entraînement du modèle pour {taxo_target} ({taxo_level}) : {e}")

    phylo_path = f"{exp_dir}/phylo_tree.txt"
    os.makedirs(os.path.dirname(phylo_path), exist_ok=True)

    with open(phylo_path, 'wb') as jtree:
        pickle.dump(phylo_tree, jtree)

    model_time = round((time.time() - start_model))
    logger.info(f"Finished make_model in {model_time} s  tree @ {phylo_path} ")



def validate(input_files, exp_dir,  params, logger,  num_processes=4, save_raw_pred=False):
    start_validation = time.time()
    logger.info(f"Start evaluation for {len(input_files)} genome files")

    metrics = ConfusionMatrixTracker()

    val_dir = os.path.join(exp_dir, 'eval')
    os.makedirs(f"{val_dir}", exist_ok=True)

    phylo_path = f"{exp_dir}/phylo_tree.txt"
    with open(phylo_path, 'rb') as jtree:
        phylo_tree: Tree = pickle.load(jtree)

    model_dir = f"{exp_dir}/model"
    process_genome_partial = partial(process_genome, phylo_tree=phylo_tree, model_dir=model_dir,
                                     params=params, val_dir=val_dir, logger=logger, metrics=metrics,
                                     num_processes=num_processes, save_raw_pred=save_raw_pred)


    # Parallelize across genomes
    with ThreadPoolExecutor(max_workers=num_processes) as executor:
        # Submit the processing of each genome as a task to the executor
        futures = [executor.submit(process_genome_partial, genome) for genome in input_files]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Predicting Genomes"):
            future.result()  #Handle exceptions by raising them if any

    log_val_metrcis(metrics, val_dir, logger)
    validation_time = round((time.time() - start_validation))
    logger.info(f"Finished validation  in {validation_time} s")


def process_genome(genome, phylo_tree, model_dir, params, val_dir, logger, metrics, num_processes, save_raw_pred):
    base_name = os.path.basename(genome).split('.')[0]
    taxons = base_name.split('_')
    logger.debug('-' * 60)
    logger.debug(f" -> Predicting for file {base_name}")

    # Ensure there are at least 6 taxonomic levels
    if len(taxons) < len(TAXO_LEVELS):
        raise ValueError(f"Genome file '{genome}' must contain at least {len(TAXO_LEVELS)} "
                         f"taxonomic levels separated by underscores.")

    gt_taxons = dict(zip(TAXO_LEVELS, taxons[:6]))

    with open(genome, 'r', encoding='utf-8') as freader:
        genome_data = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')}

    sequences = [(id_sequence, dna_sequence) for id_sequence, dna_sequence in genome_data.items()]
    partial_pred = partial(prediction, tree=phylo_tree, model_dir=model_dir, params=params, val_dir=val_dir)

    prediction_results = []

    for seq_id, seq_data in sequences:
        try:
            result =  partial_pred(seq_id, seq_data)
            pred_taxons = extract_majority_classification(result)
            metrics.update(true_labels=gt_taxons, pred_labels=pred_taxons)
            logger.debug(f"seq_id {seq_id} -> pred: {pred_taxons} -> gt: {gt_taxons}")
            prediction_results.append(result)
        except Exception as e:
            logger.debug(f"⚠️ Error for id {seq_id}: {e}")
            prediction_results.append(None)

    if save_raw_pred:
        genome_name = os.path.basename(genome).rsplit('.', 1)[0]
        report_path = os.path.join(val_dir, 'raw_pred', f"{genome_name}_job_output.json")
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        for_report = {seq_id: result for (seq_id, _), result in zip(sequences, prediction_results)}
        with open(report_path, 'w', encoding='utf-8') as jwriter:
            json.dump(for_report, jwriter)


def log_val_metrcis(metrics, val_dir, logger):
    all_val_conf_matrix = metrics.get_all_confusion_matrices()
    logger.info("=" * 60)
    logger.info("VALIDATION metrics")
    for id_level, level in enumerate(TAXO_LEVELS[:-1]):
        conf_mat_level = all_val_conf_matrix[level]
        accuracy_level = compute_accuracy_from_conf_matrix_df(conf_mat_level)
        logger.info("-" * 20)
        logger.info(f" level {level}, accuracy {accuracy_level:03f}")
        if id_level < 2:
            logger.info("\n" + conf_mat_level.to_markdown())

        file_csv = os.path.join(os.path.join(val_dir, 'metrics'), f"ConfMat_{level}.csv")
        os.makedirs(os.path.dirname(file_csv), exist_ok=True)

        conf_mat_level.to_csv(file_csv, sep=';', index=True)
        plot_conf_mat(conf_mat_level, title=f"Conf Mat for level {level}", filename=file_csv.replace(".csv", ""))


if __name__=='__main__':
    logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

    params_file = "params.yaml"
    datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"
    exp_dir = os.path.abspath('../../../exp/model_laptop')
    #
    input_files = [os.path.abspath(os.path.join(dirpath, f))
                   for dirpath, _, filenames in os.walk(datadir)
                   for f in filenames]

    with open(params_file, 'r') as file:
        params = yaml.safe_load(file)

    check_parameters(params)
    all_val_conf_matrix = validate(input_files, exp_dir, params, num_processes=4)

    print("Validation Metrics")
    for level in TAXO_LEVELS[:-1]:
        conf_mat_level = all_val_conf_matrix[level]
        accuracy_level = compute_accuracy_from_conf_matrix_df(conf_mat_level)
        print(f" level {level}, accuracy {accuracy_level}" )
        print(conf_mat_level.to_markdown())

