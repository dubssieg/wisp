import argparse
import os
import sys
import mlflow
import yaml
import logging
from datetime import datetime

from utils import setup_logger
from create_database import check_parameters
from training_functions import  validate

sys.path.append('../../..')
from wisp.wisp_light.dataset.bactero_set import BacteriaDataset

parser = argparse.ArgumentParser(description="Script d'entraînement pour le modèle bactérien.")
parser.add_argument("--exp_name", type=str, default="model_base", help="Nom de l'expérience.")
parser.add_argument("--datadir", type=str, default="/projects/microtaxo/data/refseq_with_taxo_merged", help="Répertoire des données.")
parser.add_argument("--params_file", type=str, default="params.yaml", help="Chemin du fichier de paramètres.")
parser.add_argument("--exp_rootdir", type=str, default=os.path.abspath('../../exp/'), help="Répertoire racine des expériences.")
parser.add_argument("--db_json", type=str, default=None, help="Fichier JSON de la base de données existante.")

args = parser.parse_args()
#args.datadir ="/projects/microtaxo/data/refseq_with_taxo_merged"
args.datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"
# args.datadir = "/home/hcourtei/Projects/MicroTaxo/codes/genouest_data/refseq_with_taxo_merged"
# args.datadir = "/home/hcourtei/Projects/MicroTaxo/codes/genouest_data/refseq_with_taxo_merged"
day_month_min = datetime.now().strftime('%m_%d_%H_%M')

exp_dir = "/home/hcourtei/Projects/MicroTaxo/codes/wisp/exp/model_base_03_07_16_10/"
#exp_dir = "/projects/microtaxo/exp_refseq/model_base_02_24_15_18"

os.makedirs(exp_dir, exist_ok=True)
log_file = f"{exp_dir}/eval_{day_month_min}.log"
mlflow.set_tracking_uri(f"sqlite:///{os.path.dirname(exp_dir)}/mlflow.db")
mlflow.set_experiment(args.exp_name)

logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=log_file)


with open(args.params_file, 'r') as file:
    params = yaml.safe_load(file)
    check_parameters(params)

params.update({'exp_name': args.exp_name,"day_month_min": day_month_min,
               'exp_dir':exp_dir,'datadir':args.datadir })

params_copy_path = os.path.join(exp_dir, "params.yaml")
with open(params_copy_path, 'w') as f:
    yaml.safe_dump(params, f)

with mlflow.start_run():
    mlflow.log_params(params)

    logger.info(f"Fichier {args.params_file} copié dans {params_copy_path}")

    dataset = BacteriaDataset(args.datadir, logger)
    dataset.filter_family_by_min_species(min_family_threshold=params['min_family_threshold'],
                                         max_family_repr=params['max_family_repr'])

    train_files_list, val_files_list = dataset.train_test_split(test_size=params['test_size'],
                                                                random_state=params['random_state'])


    # train_model_targets(phylo_tree,  exp_dir, params, logger, num_processes=params['num_processes'])
    validate(val_files_list, exp_dir, params, logger, max_workers=params['num_processes'])
