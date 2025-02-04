import os
import sys
import yaml
import logging
from datetime import datetime
from create_database import check_parameters
from utils import setup_logger
from training_functions import train, validate
sys.path.append('../../..')
from wisp.wisp_light.dataset.bactero_set import BacteriaDataset, TAXO_LEVELS
from wisp.wisp_light.training.metrics import  compute_accuracy_from_conf_matrix_df

day_month = datetime.now().strftime('%d_%m')
# datadir = "/home/genouest/cnrs_umr6074/hcourtei/micro_project/data/refseq_with_taxo_merged"
datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"
params_file = "params.yaml"
exp_rootdir =  os.path.abspath('../../exp/')
exp_name = 'model0'
exp_dir = f"{exp_rootdir}/{exp_name}_{day_month}"
os.makedirs(exp_dir, exist_ok=True)
log_file = f"{exp_name}_{day_month}.log"
logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=log_file)


with open(params_file, 'r') as file:
    params = yaml.safe_load(file)
    check_parameters(params)

dataset = BacteriaDataset(datadir)
dataset.filter_family_by_min_species(min_family_threshold=params['min_family_threshold'],
                                     max_family_repr=params['max_family_repr'])
# print(*file_names, sep='\n')
train_files_list, val_files_list = dataset.train_test_split(test_size=params['test_size'],
                                                            random_state=params['random_state'])


    
train(train_files_list, exp_dir, params, logger, num_processes=params['num_processes'])

all_val_conf_matrix = validate(val_files_list, exp_dir, params,logger, num_processes=4)


logger.info("="*60)
logger.info("VALIDATION metrics")
for level in TAXO_LEVELS[:-1]:
    conf_mat_level = all_val_conf_matrix[level]
    accuracy_level = compute_accuracy_from_conf_matrix_df(conf_mat_level)
    logger.info("-" * 20)
    logger.info(f" level {level}, accuracy {accuracy_level}" )
    logger.info("\n" + conf_mat_level.to_markdown())
