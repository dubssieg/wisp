import os
import time
import re
import yaml
from tqdm import tqdm
from json import load, dump
from pickle import dump as pdump, load as pload
from concurrent.futures import ThreadPoolExecutor
from create_database import build_database, validate_parameters
from create_model import make_model
from functools import partial
# from tharospytools.multithreading import futures_collector
# from utils import display_json_preview

database_name = "refseq"
input_folder = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo"
params_file = "./params.yaml"
output_dir = os.path.abspath('../../output_dir')
num_processes = 4

def is_base_file(filename):
    # Vérifie si le nom du fichier se termine par ".fna" sans numéro avant l'extension
    # et qu'il contient exactement 6 niveaux de classif séparés par des underscores
    return bool(re.match(r"^([a-zA-Z]+_){5}[a-zA-Z]+\.fna$", filename))

with open(params_file, 'r') as file:
    params = yaml.safe_load(file)

if not validate_parameters(params['DATA']):
    raise RuntimeError("Incorrect parameter file")

print("Starting database creation")
start_database = time.time()
all_paths = [os.path.abspath(os.path.join(dirpath, f))
                    for dirpath, _, filenames in os.walk(input_folder)
                        for f in filenames if is_base_file(f)
]
database_json = f'{output_dir}/databases/{database_name}.json'
os.makedirs(os.path.dirname(database_json), exist_ok=True)
phylo_tree = build_database(database_json, params['DATA'], database_name, all_paths)

database_time = round((time.time() - start_database) / 60)

print(f"Database successfully built @ {database_json} in {database_time} min")
start_model = time.time()
levels = ['root', 'domain', 'phylum', 'group', 'order']
nodes_per_level: dict = {level: [node.tag for node in list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))]
                         for i, level in enumerate(levels)} # jusqu'à order (drop family level)

# display_json_preview(output_file, num_elements=1)
print("Starting model creation")
# Loading data => should be put in the main call to escape loading it at each iteration
with open(database_json, 'r', encoding='utf-8') as jdb:
    datas = load(jdb)

# datas = {'datas': list_59_data,'mappings': taxa_code_by_level}
classif_targets = [(taxo_level, taxo_target)
                     for taxo_level, targets in nodes_per_level.items()
                     for taxo_target in targets
                     ]

make_model_partial = partial(make_model, output_dir, datas, database_json, params['MODEL'])

with ThreadPoolExecutor(max_workers=num_processes) as executor:
    futures = [executor.submit(make_model_partial, *classif_target) for classif_target in classif_targets]

for future in futures:  # Afficher une barre de progression
    model_path, config_path = future.result()  # Récupérer les résultats du future
    task_index = futures.index(future)  # accès aux arguments correspondant
    taxonomic_level, target_taxa = classif_targets[task_index]

    if model_path is not None and config_path is not None:
        key = f"{target_taxa.lower()}_{taxonomic_level}"
        try:
            node = phylo_tree[key]
            node.data.model_path = os.path.basename(model_path)  # CHANGE to have relative path
            node.data.config_path = os.path.basename(config_path)

        except KeyError:
            print(f"KEY error for phylo tree  key {key} ")
            phylo_tree.remove_node(target_taxa.lower())

phylo_path = f"{output_dir}/model/{database_name}_phylo_tree.txt"
os.makedirs(os.path.dirname(phylo_path), exist_ok=True)

with open(phylo_path, 'wb') as jtree:
    pdump(phylo_tree, jtree)

model_time = round((time.time() - start_model) / 60)
print(f"range]Finished computing models, tree @ {phylo_path} in {model_time} min : TOTAL {database_time + model_time} min")





