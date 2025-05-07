import logging
import os
import time
import sys
import pickle
import json
from tqdm import tqdm
from treelib import Tree
from Bio import SeqIO
import mlflow

from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from wisp_light.training.create_model import make_model
from wisp_light.training.create_prediction import prediction
from wisp_light.training.utils import  extract_majority_classification, setup_logger
from wisp_light.training.metrics import ConfusionMatrixTracker, compute_accuracy_from_conf_matrix_df

from wisp_light.visu.plots_tools import plot_conf_mat
from wisp_light.dataset.refSeqDataset import TAXO_LEVELS
from wisp_light.training.utils import log_resource_usage



def train_model_targets(phylo_tree, exp_dir, params, logger, max_workers=4):
    """
    Entraîne un modèle pour chaque nœud du niveau taxonomique donné dans un arbre phylogénétique.

    Parameters
    ----------
    phylo_tree : Tree
        Arbre phylogénétique Treelib avec les taxons et leurs relations hiérarchiques.
    exp_dir : str
        Répertoire de sauvegarde des modèles et du fichier phylo_tree sérialisé.
    params : dict
        Dictionnaire de paramètres d'entraînement (extraits d'un fichier YAML).
    logger : logging.Logger
        Logger pour affichage des messages d'information, debug et erreurs.
    max_workers : int, optional
        Nombre maximal de threads pour l'entraînement parallèle, par défaut 4.

    Returns
    -------
    int
        Temps total d'entraînement en secondes.
    """

    start_model = time.time()

    nodes_per_level: dict = {
        level: [node.tag for node in list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))]
        for i, level in enumerate(['root']+ TAXO_LEVELS[:-1])}  # jusqu'à order (drop family level)

    # display_json_preview(output_file, num_elements=1)
    classif_targets = [(taxo_level, taxo_target)
                       for taxo_level, targets in nodes_per_level.items()
                       for taxo_target in targets
                       ]

    with open(f'{exp_dir}/databases.json', 'r', encoding='utf-8') as jdb:
        database = json.load(jdb)  # Loading data => should be put in the main call to escape loading it at each iteration

    mappings_data = database['mappings']
    all_data = database['datas']

    logger.info("Starting model creation")

    def get_filtered_reads(all_data, taxo_level, taxo_target):
        return [
            (data_by_genome, read)
            for data_by_genome in all_data
            if taxo_level == "root" or data_by_genome.get(taxo_level) == taxo_target
            for read in data_by_genome['datas']
        ]
    # make_model_partial = partial(make_model, exp_dir, database, params, logger)
    logger.info(f"Lancement de {len(classif_targets)} modèles avec num_processes={max_workers}")
    nb_model_fail = 0  # Compteur de modèles non générés

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                make_model,
                output_dir=exp_dir,
                filtered_reads=get_filtered_reads(all_data, taxo_level, taxo_target),
                mappings_data=mappings_data,
                params=params,
                logger=logger,
                taxo_level=taxo_level,
                taxo_target=taxo_target
            ): (taxo_level, taxo_target)
            for taxo_level, taxo_target in classif_targets
        }
        log_resource_usage(logger, "train_model_targets (Pool)")

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
                    nb_model_fail += 1
            except Exception as e:
                logger.error(f"🔥 Erreur lors de l'entraînement du modèle pour {taxo_target} ({taxo_level}) : {e}")

        logger.info(f"📊 Modèles non générés : {nb_model_fail}/{ len(futures) }")

    phylo_path = f"{exp_dir}/phylo_tree.txt"
    os.makedirs(os.path.dirname(phylo_path), exist_ok=True)

    with open(phylo_path, 'wb') as jtree:
        pickle.dump(phylo_tree, jtree)

    model_time = round((time.time() - start_model))
    logger.info(f"Finished make_model in {model_time} s  tree @ {phylo_path} ")
    mlflow.log_metric("model_time", model_time)
    log_resource_usage(logger, "train_model_targets (end)")
    return model_time



def validate(val_dataset, exp_dir, params, logger, max_workers=4, save_raw_pred=False):
    """
    Valide les modèles entraînés sur un jeu de données de validation.

    Parameters
    ----------
    val_dataset : list of tuple
        Liste de tuples (path_to_genome_file, dict_gt_taxons), un par génome à valider.
    exp_dir : str
        Répertoire de l'expérience contenant les modèles et le phylo_tree.
    params : dict
        Paramètres utilisés pour la prédiction.
    logger : logging.Logger
        Logger pour le suivi de l'exécution.
    max_workers : int, optional
        Nombre de processus parallèles pour la prédiction, par défaut 4.
    save_raw_pred : bool, optional
        Si True, sauvegarde les prédictions brutes au format JSON, par défaut False.

    Returns
    -------
    int
        Durée de la validation en secondes.
    """

    start_validation = time.time()
    logger.info(f"Start evaluation for {len(val_dataset)} genome files")

    metrics = ConfusionMatrixTracker()

    val_dir = os.path.join(exp_dir, 'eval')
    os.makedirs(f"{val_dir}", exist_ok=True)

    phylo_path = f"{exp_dir}/phylo_tree.txt"
    with open(phylo_path, 'rb') as jtree:
        phylo_tree: Tree = pickle.load(jtree)

    model_dir = os.path.join(exp_dir,"model")
    process_genome_partial = partial(process_genome, phylo_tree=phylo_tree, model_dir=model_dir,
                                     params=params, val_dir=val_dir, logger=logger,
                                     save_raw_pred=save_raw_pred)

    all_metrics_samples = []  # Liste pour stocker les métriques de chaque échantillon
    # Parallelize across genomes
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit the processing of each genome as a task to the executor
        futures = [executor.submit(process_genome_partial, sample) for sample in val_dataset]
        log_resource_usage(logger, "validate (Pool)")

        for future in tqdm(as_completed(futures), total=len(futures), desc="Predicting Genomes", disable=not sys.stdout.isatty()):
            metrics_sample = future.result()  #Handle exceptions by raising them if any
            all_metrics_samples.extend(metrics_sample)

    for true_labels, pred_labels in all_metrics_samples:
        metrics.update(true_labels=true_labels, pred_labels=pred_labels)

    log_val_metrics(metrics, val_dir, logger)
    validation_time = round((time.time() - start_validation))
    logger.info(f"Finished validation  in {validation_time} s")
    mlflow.log_metric("validation_time", validation_time)
    log_resource_usage(logger, "validate (end)")
    return validation_time


def process_genome(sample, phylo_tree, model_dir, params, val_dir, logger, save_raw_pred):
    """
    Prédit les taxons d'un génome à l'aide de l'arbre phylogénétique et des modèles associés.

    Parameters
    ----------
    sample : tuple
        Tuple (genome_path, ground_truth_taxons) pour un échantillon.
    phylo_tree : Tree
        Arbre phylogénétique contenant les chemins vers les modèles.
    model_dir : str
        Répertoire contenant les modèles entraînés.
    params : dict
        Paramètres utilisés pour la prédiction.
    val_dir : str
        Répertoire dans lequel sont stockés les résultats de validation.
    logger : logging.Logger
        Logger pour le suivi des événements et erreurs.
    save_raw_pred : bool
        Si True, sauvegarde les prédictions individuelles.

    Returns
    -------
    list of tuple
        Liste de tuples (true_labels, predicted_labels) pour chaque séquence du génome.
    """

    genome, gt_taxons = sample
    metrics_sample = []
    # base_name = os.path.basename(genome).split('.')[0]
    # taxons = base_name.split('_')
    logger.debug('-' * 60)
    logger.debug(f" -> Predicting for file {genome}")

    # Ensure there are at least 6 taxonomic levels
    if not list(gt_taxons.keys()) == TAXO_LEVELS :
        raise ValueError(f"Genome file '{genome}' must contain at least {len(TAXO_LEVELS)} "
                         f"taxonomic levels separated by underscores.")

    with open(genome, 'r', encoding='utf-8') as freader:
        genome_data = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')}

    logger.debug(f"File {os.path.basename(genome)} contains {len(genome_data)} sequences.")
    sequences = [(id_sequence, dna_sequence) for id_sequence, dna_sequence in genome_data.items()]
    partial_pred = partial(prediction, tree=phylo_tree, model_dir=model_dir, params=params, val_dir=val_dir, logger=logger)

    prediction_results = []
    file_error_plylum = open(f"{val_dir}/error_phylum_val.txt", "a")
    for seq_id, seq_data in sequences:
        try:
            result =  partial_pred(seq_id, seq_data)
            logger.debug(f"[DEBUG] Prediction result for {seq_id}: {result}")
            pred_taxons = extract_majority_classification(result)
            if not gt_taxons:
                logger.warning(f"⚠️ Empty ground-truth taxons for {genome}")
            metrics_sample.append((gt_taxons, pred_taxons))

            if gt_taxons['phylum'] != pred_taxons['phylum']:
                logger.debug("ERROR phylum")
                logger.debug(f"-> seq_id {seq_id} -> pred: {pred_taxons} -> gt: {gt_taxons}")
                file_error_plylum.write(f"seq_id {seq_id}\n pred: {pred_taxons}\n gt  : {gt_taxons} \n")

            prediction_results.append(result)
        except Exception as e:
            if not isinstance(e, ValueError):
                logger.exception(f"⚠️ Exception for sequence ID {seq_id}: {e}")
            prediction_results.append(None)
    file_error_plylum.close()
    if save_raw_pred:
        genome_name = os.path.basename(genome).rsplit('.', 1)[0]
        report_path = os.path.join(val_dir, 'raw_pred', f"{genome_name}_job_output.json")
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        for_report = {seq_id: result for (seq_id, _), result in zip(sequences, prediction_results)}
        with open(report_path, 'w', encoding='utf-8') as jwriter:
            json.dump(for_report, jwriter)

    return metrics_sample


def log_val_metrics(metrics, val_dir, logger):
    """
    Calcule et enregistre les métriques de validation pour chaque niveau taxonomique.

    Parameters
    ----------
    metrics : ConfusionMatrixTracker
        Objet de suivi des matrices de confusion et des métriques par niveau.
    val_dir : str
        Répertoire où sauvegarder les métriques et graphiques.
    logger : logging.Logger
        Logger pour afficher les résultats.
    """
    # all_val_conf_matrix = metrics.get_all_confusion_matrices()
    logger.info("=" * 60)
    logger.info("VALIDATION metrics")
    # print("DEBUG - Keys in true_labels:", metrics.true_labels.keys())
    # for level in TAXO_LEVELS:
    #     print(f"Level: {level}, Nb labels: {len(metrics.true_labels[level])}")
    metrics.build_taxonomy_df()
    LEVEL_MARKER = 'phylum'

    for id_level, level in enumerate(TAXO_LEVELS):
        conf_mat_level = metrics.get_confusion_matrix(level)
        accuracy_level = compute_accuracy_from_conf_matrix_df(conf_mat_level)
        mlflow.log_metric(f"accuracy_{level}", accuracy_level)
        logger.info("-" * 20)
        logger.info(f" level {level}, accuracy {accuracy_level:03f}")

        if len(conf_mat_level) <= 15 :
            logger.info("\n" + conf_mat_level.to_markdown())

        file_csv = os.path.join(os.path.join(val_dir, 'metrics'), f"ConfMat_{level}.csv")
        os.makedirs(os.path.dirname(file_csv), exist_ok=True)

        conf_mat_level.to_csv(file_csv, sep=';', index=True)
        mlflow.log_artifact(file_csv)

        plot_path = file_csv.replace(".csv", ".png")
        if level != LEVEL_MARKER:
            separator_indices = metrics.calculate_separator_indices(level_marker=LEVEL_MARKER, level_index=level)
        else:
            separator_indices = None

        plot_conf_mat(conf_mat_level, level, separator_indices, filename=plot_path)
        mlflow.log_artifact(plot_path)

def count_seq(dataset):
    """
    Compte le nombre total de séquences dans un jeu de données.

    Parameters
    ----------
    dataset : list of tuple
        Liste de tuples (genome_path, ground_truth_taxons).

    Returns
    -------
    int
        Nombre total de séquences.
    """
    nb_seq = 0
    for sample in dataset:
        genome, gt_taxons = sample
        with open(genome, 'r', encoding='utf-8') as freader:
            genome_data = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')}

        sequences = [(id_sequence, dna_sequence) for id_sequence, dna_sequence in genome_data.items()]
        for seq_id, seq_data in sequences:
            nb_seq += 1

    print(f" nb sequence {nb_seq} for nb val genome {len(dataset)}")
    return nb_seq

if __name__=='__main__':
    logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

    params_file = "params.yaml"
    datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"
    exp_dir = os.path.abspath('../../../exp/model_laptop')
