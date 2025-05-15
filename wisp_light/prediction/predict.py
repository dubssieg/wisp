import os
import pickle
import logging
import yaml
from pprint import pformat
from Bio import SeqIO

from wisp_light.training.create_prediction import prediction
from wisp_light.training.utils import setup_logger, extract_majority_classification

log_file = None
exp_dir = "/home/hcourtei/Projects/MicroTaxo/codes/exp/model_base_complete_05_15_16_37"
save_raw_pred = False
val_dir = "/home/hcourtei/Projects/MicroTaxo/codes/predict"

logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=os.path.join(val_dir,'logs'))

phylo_path = f"{exp_dir}/phylo_tree.txt"
with open(phylo_path, 'rb') as jtree:
    phylo_tree = pickle.load(jtree)

params_file = "/home/hcourtei/Projects/MicroTaxo/codes/wisp_light/training/params.yaml"

with open(params_file, 'r') as file:
    params = yaml.safe_load(file)

model_dir = os.path.join(exp_dir, "model")


genome = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data/GCF_000725405.1_ASM72540v1_genomic.fna"


with open(genome, 'r', encoding='utf-8') as freader:
    genome_data = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')}

print(f"File {os.path.basename(genome)} contains {len(genome_data)} sequences.")
sequences = [(id_sequence, dna_sequence) for id_sequence, dna_sequence in genome_data.items()]
sequence = sequences[0]
id_sequence, dna_sequence = sequence

results = prediction(id_sequence, dna_sequence, params, phylo_tree, model_dir, val_dir, logger)
pred_taxons = extract_majority_classification(results)

logger.info("="*60)
logger.info(f"Raw prediction for genome 0 of {os.path.basename(genome)}")
logger.info(f"with xgboost model from model dir {exp_dir}")
logger.info('\n'+pformat(results))
logger.info("-"*60)
logger.info("Majority prediction")
logger.info('\n'+pformat(pred_taxons))
logger.info("="*60)
