"""
Script d'entraînement pour le modèle de classification bactérienne.

Ce script permet de :
- Créer ou charger une base de données phylogénétique à partir d’un ensemble de génomes.
- Entraîner un modèle de classification basé sur cette base.
- Valider le modèle sur un jeu de validation.
- Suivre et journaliser les métriques avec MLflow.
- Gérer les logs et la configuration via des fichiers YAML.

Utilisation :
-------------
Ce script s'utilise en ligne de commande avec plusieurs arguments optionnels :
    --exp_name : nom de l'expérience (défaut = "model_base_complete")
    --datadir : répertoire contenant les génomes FASTA (défaut = refseq_with_taxo_merged)
    --params_file : chemin vers le fichier YAML des hyperparamètres (défaut = params.yaml)
    --db_json : chemin vers une base de données préexistante à recharger (optionnel)

Exemple :
---------
python train.py --exp_name test_05 --params_file config.yaml --db_json path/to/db.json

Modules utilisés :
------------------
- Bio.SeqIO : pour parser les fichiers FASTA
- argparse : gestion des arguments en ligne de commande
- yaml : chargement des paramètres
- logging : journalisation des événements
- mlflow : suivi des expériences
- wisp_light.dataset.refSeqDataset : gestion du dataset bactérien
- training_functions : fonctions d’entraînement et de validation
- create_database : construction et chargement de la base phylogénétique

Auteur : Hermann Courteille
"""
import argparse
import os
import time
import yaml
import logging
import mlflow
from pathlib import Path
from datetime import datetime
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

from wisp_light.training.utils import setup_logger
from wisp_light.training.create_database import  build_database, load_phylo_tree
from wisp_light.training.training_functions import train_model_targets, validate, count_seq

from wisp_light.dataset.refSeqDataset import RefSeqDataset
# from wisp.wisp_light.dataset.bactero_set import BacteriaDataset

parser = argparse.ArgumentParser(description="Script d'entraînement pour le modèle bactérien.")
parser.add_argument("--exp_name", type=str, default="model_base_complete", help="Nom de l'expérience.")
parser.add_argument("--datadir", type=str, default="/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data", help="Répertoire des données.")
parser.add_argument("--params_file", type=str, default="training/params.yaml", help="Chemin du fichier de paramètres.")
parser.add_argument("--index_csv", type=str, default="build_dataset/refseq/reference_genome_summary_complete_taxo.tsv", help="Chemin vers le fichier CSV d'index des génomes.")
parser.add_argument("--db_json", type=str, default="", help="Fichier JSON de la base de données déjà existante.")



args = parser.parse_args()
wisplight_dir =  Path(__file__).resolve().parent.parent
print("ROOT wisp_light", wisplight_dir)
index_csv = wisplight_dir/ args.index_csv
params_file = wisplight_dir / args.params_file
exp_rootdir =  wisplight_dir.parent / 'exp'

# args.db_json = '/home/hcourtei/Projects/MicroTaxo/codes/exp_refseq/model_base_index_03_24_15_06/databases.json'

# args.datadir = '/WORKS/microtaxo/data/refseq3' # '/scratch/hcourtei/refseq3'
# args.exp_rootdir = '/WORKS/microtaxo/exp_refseq' # '/scratch/hcourtei/exp_refseq'
# args.db_json  = "/WORKS/microtaxo/exp_refseq/model_base_testval_04_04_13_48/databases.json"


day_month_min = datetime.now().strftime('%m_%d_%H_%M')
if args.db_json:
    exp_dir = os.path.dirname(args.db_json)
    log_file = f"{exp_dir}/restart.log"

else:

    exp_dir = f"{exp_rootdir}/{args.exp_name}_{day_month_min}"
    os.makedirs(exp_dir, exist_ok=True)
    log_file = f"{exp_dir}/init_train.log"


logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=log_file)

mlflow.set_tracking_uri(f"file://{os.path.dirname(exp_dir)}/mlruns") # "file://chemin_ml_runs"
logger.info(f"Current tracking uri: { mlflow.get_tracking_uri()}")
mlflow.set_experiment(args.exp_name)

logger.info(f"current working directory: {os.getcwd()}")

with open(params_file, 'r') as file:
    params = yaml.safe_load(file)

params.update({'exp_name': args.exp_name,"day_month_min": day_month_min,
               'exp_dir':exp_dir,'datadir':args.datadir, 'db_json': args.db_json})

params_copy_path = os.path.join(exp_dir, "params.yaml")
with open(params_copy_path, 'w') as f:
    yaml.safe_dump(params, f)


print(f" nb core cpu {os.cpu_count()} , counting kmer with max_workers_trainval {params['max_workers_trainval']} "
      f"max_workers_db {params['max_workers_db']}")

with mlflow.start_run():
    mlflow.log_params(params)

    logger.info(f"Fichier {params_file} copié dans {params_copy_path}")

    # dataset = BacteriaDataset(args.datadir, logger)
    # dataset.filter_family_by_min_species(min_family_threshold=params['min_family_threshold'],
    #                                      max_family_repr=params['max_family_repr'])
    dataset = RefSeqDataset(index_csv, args.datadir, logger, cut= params['cut_data'])

    train_dataset, val_dataset = dataset.split(test_size=params['test_size'], random_state=params['random_state'],
                                               family_strat=params['family_strat'])

    nb_seq = count_seq(val_dataset)

if  args.db_json:

    phylo_tree, nb_genome_indexed = load_phylo_tree(args.db_json)
    logger.info(f"Reload database json {nb_genome_indexed} genomes indexed from {exp_dir} ")

else:
    logger.info(f"Starting database creation for {len(train_dataset)} genome files ")
    start_database = time.time()
    database_json = os.path.join(exp_dir, f"databases.json")
    phylo_tree = build_database(train_dataset, params, database_json, logger,max_workers=params['max_workers_db'])
    database_time = round((time.time() - start_database))
    logger.info(f"Database successfully built in {database_time} s @ {f'{exp_dir}/databases.json'} ")
    mlflow.log_metric("database_time", database_time)

model_time = train_model_targets(phylo_tree, exp_dir, params, logger, max_workers=params['max_workers_trainval'])
#
validation_time = validate(val_dataset, exp_dir, params, logger, max_workers=params['max_workers_trainval'])

print(f"Times: \n - database {database_time} s\n - model {model_time} s\n - validation {validation_time} s")
