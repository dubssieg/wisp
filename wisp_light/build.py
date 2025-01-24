import os
import time
import re
from tqdm import tqdm
from json import load, dump
from pickle import dump as pdump, load as pload
from concurrent.futures import ThreadPoolExecutor
# from tharospytools.multithreading import futures_collector
from create_database import build_database
from create_model import make_model

database_name = "refseq"
input_folder = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo"
parameters = "./params.yaml"
output_dir = os.path.abspath('../../output_dir')
num_processes = 4

# Fonction pour filtrer les fichiers sans numéro à la fin
def is_base_file(filename):
    # Vérifie si le nom du fichier se termine par ".fna" sans numéro avant l'extension
    # et qu'il contient exactement 6 niveaux de classif séparés par des underscores
    return bool(re.match(r"^([a-zA-Z]+_){5}[a-zA-Z]+\.fna$", filename))
    # Vérifie si le nom du fichier se termine par ".fna" sans numéro avant l'extension
    # return bool(re.match(r".*[^_\d]\.fna$", filename))


print("Starting database creation")
start_database = time.time()
all_paths = [os.path.abspath(os.path.join(dirpath, f))
                    for dirpath, _, filenames in os.walk(input_folder)
                        for f in filenames if is_base_file(f)
]

output_path, phylo_tree = build_database(output_dir, parameters, database_name, all_paths)

database_time = round((time.time() - start_database) / 60)

print(f"Database sucessfully built @ {output_path} in {database_time} min")
start_model = time.time()

nodes_per_level: dict = {level: [node.tag for node in list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))]
                         for i, level in enumerate(['root', 'domain', 'phylum', 'group', 'order'], start=0)}

print("Starting model creation")

# Loading data => should be put in the main call to escape loading it at each iteration
with open(output_path, 'r', encoding='utf-8') as jdb:
    datas = load(jdb)
task_list= [(output_dir, datas, output_path, taxonomic_level, target_taxa)
                                                         for taxonomic_level, targets in nodes_per_level.items() for
                                                         target_taxa in targets]

# retcodes: list = futures_collector(make_model, task_list)
  # Remplacez ce nombre par celui que vous souhaitez utiliser

with ThreadPoolExecutor(max_workers=num_processes) as executor:
    futures = [executor.submit(make_model, *task) for task in task_list]

for future in tqdm(futures):  # Afficher une barre de progression
    model_path, config_path = future.result()  # Récupérer les résultats du future

    # Accéder à l'index du future pour trouver le task correspondant
    task_index = futures.index(future)
    _, _, _, taxonomic_level, target_taxa = task_list[task_index]

    if model_path is not None and config_path is not None:

        try:
            key = f"{target_taxa.lower()}_{taxonomic_level}"
            node = phylo_tree[key]
            node.data.model_path = os.path.basename(model_path)  # CHANGE to have relative path
            node.data.config_path = os.path.basename(config_path)

        except KeyError:
            print(f"KEY error for phylo tree  key {key} ")
            phylo_tree.remove_node(target_taxa.lower())

# phylo_path = "/home/hcourtei/Projects/MicroTaxo/codes/envwisp/lib/python3.12/site-packages/workspace/model/resfseq_phylo_tree.txt"
phylo_path = f"{output_dir}/model/{database_name}_phylo_tree.txt"
os.makedirs(os.path.dirname(phylo_path), exist_ok=True)

with open(phylo_path, 'wb') as jtree:
    pdump(phylo_tree, jtree)

model_time = round((time.time() - start_model) / 60)
print(f"range]Finished computing models, tree @ {phylo_path} in {model_time} min : TOTAL {database_time + model_time} min")
