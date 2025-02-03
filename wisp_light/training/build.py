import os
import time
import yaml
import logging
from json import load
from pickle import dump as pdump
from concurrent.futures import ThreadPoolExecutor
from create_database import build_database, check_parameters
from create_model import make_model
from functools import partial
from utils import setup_logger
from wisp.wisp_light.dataset.bactero_set import BacteriaDataset



logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

database_name = "refseq"
params_file = "training/params.yaml"
output_dir = os.path.abspath('../../exp/model_laptop')
num_processes = 4
datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"
print(f"outputdir {output_dir}")



with open(params_file, 'r') as file:
    params = yaml.safe_load(file)


check_parameters(params)

dataset = BacteriaDataset(datadir)
dataset.filter_family_by_min_species(min_family_threshold=params['min_family_threshold'],
                                     max_family_repr=params['max_family_repr'])
# print(*file_names, sep='\n')
train_files_list, val_files_list = dataset.train_test_split(test_size=params['test_size'],
                                                            random_state=params['random_state'])


logger.info("Starting database creation")
start_database = time.time()

database_json = f'{output_dir}/databases/{database_name}.json'
os.makedirs(os.path.dirname(database_json), exist_ok=True)



phylo_tree = build_database(train_files_list, params, database_json)

database_time = round((time.time() - start_database) / 60)

logger.info(f"Database successfully built @ {database_json} in {database_time} min")
start_model = time.time()

levels = ['root', 'domain', 'phylum', 'group', 'order']
nodes_per_level: dict = {level: [node.tag for node in list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))]
                         for i, level in enumerate(levels)} # jusqu'à order (drop family level)

# display_json_preview(output_file, num_elements=1)
logger.info("Starting model creation")
with open(database_json, 'r', encoding='utf-8') as jdb:
    datas = load(jdb) # Loading data => should be put in the main call to escape loading it at each iteration

# datas = {'datas': list_59_data,'mappings': taxa_code_by_level}
classif_targets = [(taxo_level, taxo_target)
                     for taxo_level, targets in nodes_per_level.items()
                     for taxo_target in targets
                     ]

make_model_partial = partial(make_model, output_dir, datas, database_json, params)

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

phylo_path = f"{output_dir}/model/{database_name}_phylo_tree.txt"
os.makedirs(os.path.dirname(phylo_path), exist_ok=True)

with open(phylo_path, 'wb') as jtree:
    pdump(phylo_tree, jtree)

model_time = round((time.time() - start_model) / 60)
logger.info(f"range]Finished computing models, tree @ {phylo_path} in {model_time} min : "
            f"TOTAL {database_time + model_time} min")






