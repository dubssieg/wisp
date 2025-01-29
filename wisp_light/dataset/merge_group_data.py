import os
import re
import sys
import shutil
import logging
import argparse
from collections import defaultdict
from tqdm import tqdm
from bactero_set import TAXO_LEVELS
sys.path.append('../..')
from wisp.wisp_light.utils import setup_logger

parser = argparse.ArgumentParser(description="merge des répertoire groupes pour refseq")
parser.add_argument("--datadir", default="/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo",)
parser.add_argument("--debug", action="store_true", help="Activer le mode debug (niveau de log DEBUG)")
args = parser.parse_args()

level = logging.DEBUG if args.debug else logging.INFO
logger = setup_logger(os.path.basename(__file__), level=level)

output_dir = f"{args.datadir}_merged" # Répertoire de sortie
os.makedirs(output_dir, exist_ok=True)

pattern_parts = [f"(?P<{level}>\\w+?)" for level in TAXO_LEVELS]
pattern_filename = "^" + "_".join(pattern_parts) + "(?:_(?P<id>\\d+))?\\.fna$"


file_counter = defaultdict(int) # Dictionnaire pour suivre les occurrences des noms de base

taxo_levels_unmatched = []
file_no_fna = []
for root, groups , files in tqdm(os.walk(args.datadir)): # Parcours des fichiers dans tous les sous-répertoires

    for file in files:
        base_name, ext = os.path.splitext(file)
        src_path = os.path.join(root, file)

        if ext == ".fna":
            match = re.match(pattern_filename, file)
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
                logger.debug(f"Copié : {os.path.basename(src_path)} → {os.path.basename(dst_path)}")
                shutil.copy2(src_path, dst_path)

            else:
                logger.debug(f"no match with {'_'.join(TAXO_LEVELS)}*(id).fna")
                taxo_levels_unmatched.append(src_path)

        else :  #if ext == ".fna":
            file_no_fna.append(src_path)

if taxo_levels_unmatched:
    with open(f"{args.datadir}_taxo_levels_unmatched.txt", "w") as f_taxo:
        f_taxo.writelines(f"{line}\n" for line in taxo_levels_unmatched)

if file_no_fna:
    with open(f"{args.datadir}/file_no_fna.txt", "w") as f_no_fna:
        f_no_fna.writelines(f"{line}\n" for line in file_no_fna)