import os
import re
import sys
import shutil
import logging
import argparse
from collections import defaultdict
from tqdm import tqdm
sys.path.append("../..")
from wisp.wisp_light.dataset.bactero_set import TAXO_LEVELS
from wisp.wisp_light.training.utils import setup_logger

parser = argparse.ArgumentParser(description="merge des répertoire groupes pour refseq")
parser.add_argument("--datadir", default="/projects/microtaxo/data/refseq_with_taxo",)
parser.add_argument("--debug", action="store_true", help="Activer le mode debug")
args = parser.parse_args()

level = logging.DEBUG if args.debug else logging.INFO
logger = setup_logger(os.path.basename(__file__), level=level, log_file=None)

output_dir = f"{args.datadir}_merged" # Répertoire de sortie
os.makedirs(output_dir, exist_ok=True)

TAXO_LEVELS = ["domain", "phylum", "group", "order", "family", "specie"]
pattern_parts = [f"(?P<{level}>\\w+?)" for level in TAXO_LEVELS]
pattern_filename = "^" + "_".join(pattern_parts) + "(?:_(?P<id>\\d+))?\\.fna$"  #A_B_C_D_E_F_***.fna

file_counter = defaultdict(int) # Dictionnaire pour suivre les occurrences des noms de base

taxo_levels_unmatched = []
for root, groups , filenames in tqdm(os.walk(args.datadir)): # Parcours des fichiers dans tous les sous-répertoires

    for filename in filenames:
        base_name, ext = os.path.splitext(filename)
        src_path = os.path.join(root, filename)
        match = re.match(pattern_filename, filename)
        if match:
            result = match.groupdict()
            base_name = "_".join(result[level] for level in TAXO_LEVELS)
            # base_name = f"{result['domain']}_{result['phylum']}_{result['group']}_{result['order']}_{result['family']}_{result['specie']}"
            if base_name in file_counter:
                file_counter[base_name] += 1
            else:
                file_counter[base_name] = 0
            new_id = file_counter[base_name]
            new_name = f"{base_name}_{new_id}{ext}"

            dst_path = os.path.join(output_dir, new_name)
            logger.debug(f"Moved and renamed : {os.path.basename(dst_path)}")
            shutil.copy2(src_path, dst_path)

        else: # if no match pattern
            logger.debug(f"WARNING no pattern match for {src_path}")
            taxo_levels_unmatched.append(src_path)

if taxo_levels_unmatched:
    with open(f"{args.datadir}_taxo_levels_unmatched.txt", "w") as f_taxo:
        f_taxo.writelines(f"{line}\n" for line in taxo_levels_unmatched)
