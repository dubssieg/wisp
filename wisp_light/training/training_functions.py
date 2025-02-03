"""Builds predictions from reads"""
import logging
import os
import yaml
import time
from tqdm import tqdm
from json import dump
from treelib import Tree
from Bio import SeqIO
from pickle import load as pload
from json import load
from pickle import dump as pdump
from create_model import make_model
from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed

from create_database import check_parameters, build_database
from create_prediction import prediction
from utils import  extract_majority_classification, setup_logger

from wisp.wisp_light.training.metrics import ConfusionMatrixTracker, compute_accuracy_from_conf_matrix_df

logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

LEVELS_gt = ['domain', 'phylum', 'group', 'order', 'family', 'species']


def train(train_files_list, exp_dir, params, num_processes=4):
    logger.info(f"Starting database creation for {len(train_files_list)} genome files ")

    start_database = time.time()

    phylo_tree = build_database(train_files_list, params, f'{exp_dir}/databases.json')

    database_time = round((time.time() - start_database) / 60)

    logger.info(f"Database successfully built @ {f'{exp_dir}/databases.json'} in {database_time} min")
    start_model = time.time()

    levels = ['root', 'domain', 'phylum', 'group', 'order']
    nodes_per_level: dict = {
        level: [node.tag for node in list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))]
        for i, level in enumerate(levels)}  # jusqu'à order (drop family level)

    # display_json_preview(output_file, num_elements=1)
    logger.info("Starting model creation")
    with open(f'{exp_dir}/databases.json', 'r', encoding='utf-8') as jdb:
        datas = load(jdb)  # Loading data => should be put in the main call to escape loading it at each iteration

    # datas = {'datas': list_59_data,'mappings': taxa_code_by_level}
    classif_targets = [(taxo_level, taxo_target)
                       for taxo_level, targets in nodes_per_level.items()
                       for taxo_target in targets
                       ]

    make_model_partial = partial(make_model, exp_dir, datas, params)

    with ThreadPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(make_model_partial, *classif_target) for classif_target in classif_targets]

    for future in futures:  # Afficher une barre de progression
        model_path = future.result()  # Récupérer les résultats du future
        task_index = futures.index(future)  # accès aux arguments correspondant
        taxonomic_level, target_taxa = classif_targets[task_index]

        if model_path is not None:
            key = f"{target_taxa.lower()}_{taxonomic_level}"
            try:
                node = phylo_tree[key]
                node.data.model_path = os.path.basename(model_path)  # CHANGE to have relative path

            except KeyError:
                logger.error(f"KEY error for phylo tree  key {key}  remove node {target_taxa.lower()}")
                phylo_tree.remove_node(target_taxa.lower())

    phylo_path = f"{exp_dir}/phylo_tree.txt"
    os.makedirs(os.path.dirname(phylo_path), exist_ok=True)

    with open(phylo_path, 'wb') as jtree:
        pdump(phylo_tree, jtree)

    model_time = round((time.time() - start_model) / 60)
    logger.info(f"range]Finished computing models, tree @ {phylo_path} in {model_time} min : "
                f"TOTAL {database_time + model_time} min")



def validate(input_files, exp_dir,  params, num_processes=4):
    start_validation = time.time()
    logger.info(f"Start evaluation for {len(input_files)} genome files")

    metrics = ConfusionMatrixTracker()

    predict_dir = os.path.join(exp_dir, 'predict')
    os.makedirs(f"{predict_dir}", exist_ok=True)

    phylo_path = f"{exp_dir}/phylo_tree.txt"
    with open(phylo_path, 'rb') as jtree:
        phylo_tree: Tree = pload(jtree)

    model_dir = f"{exp_dir}/model"


    for id_g, genome in enumerate(input_files[:10]):
        base_name = os.path.basename(genome).split('.')[0]
        taxons = base_name.split('_')
        logger.debug('-'*60)
        logger.debug(f" -> id_g {id_g} Predict file {base_name}")

        # Vérifier qu'il y a au moins 6 niveaux taxonomiques
        if len(taxons) < len(LEVELS_gt):
            raise ValueError(
                f"Le fichier '{genome}' doit contenir au moins {len(LEVELS_gt)} niveaux taxonomiques séparés par des underscores.")

        gt_taxons = dict(zip(LEVELS_gt, taxons[:6]))

        with open(genome, 'r', encoding='utf-8') as freader:
            genome_data = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')}

        sequences = [(id_sequence, dna_sequence) for id_sequence, dna_sequence in genome_data.items()]
        partial_pred = partial(prediction, tree=phylo_tree, model_dir=model_dir, params=params,output_dir=predict_dir)

        with ThreadPoolExecutor(max_workers=num_processes) as executor:
            future_to_seq = {executor.submit(partial_pred, *seq): seq for seq in sequences}
        prediction_results = []
        for future in as_completed(future_to_seq):  # Itère sur les futures terminés (pas d'ordre garanti)
            seq_id, _ = future_to_seq[future]  # Récupérer l'ID de la séquence associée

            try:
                result = future.result()
                pred_taxons = extract_majority_classification(result)
                metrics.update(true_labels=gt_taxons, pred_labels=pred_taxons)
                logger.debug(f"seq_id {seq_id}") # \n {result}")
                logger.debug(f"-> pred : {pred_taxons}")
                logger.debug(f"-> gt : {gt_taxons}")
                prediction_results.append(result)  # Récupérer le résultat si pas d'erreur
            except Exception as e:
                logger.debug(f"⚠️ Error for id {seq_id}: {e}")  # Afficher l'erreur sans arrêter
                prediction_results.append(None)  # Insérer un résultat par défaut

        genome_name = genome.split('/')[-1].rsplit('.', 1)[0]
        report_path = os.path.join(predict_dir, f"{genome_name}_job_output.json")

        for_report = {seq_id: result for (seq_id, _), result in zip(sequences, prediction_results)}
        with open(report_path, 'w', encoding='utf-8') as jwriter:
            dump(for_report, jwriter)

    all_val_conf_matrix = metrics.get_all_confusion_matrices()
    validation_time = round((time.time() - start_validation) / 60)
    logger.info(f"Finished validation  in {validation_time} min")

    return all_val_conf_matrix


if __name__=='__main__':
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
    for level in LEVELS_gt[:-1]:
        conf_mat_level = all_val_conf_matrix[level]
        accuracy_level = compute_accuracy_from_conf_matrix_df(conf_mat_level)
        print(f" level {level}, accuracy {accuracy_level}" )
        print(conf_mat_level.to_markdown())
