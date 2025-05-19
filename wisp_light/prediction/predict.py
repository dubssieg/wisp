import os
import pickle
import logging
import yaml
import pandas as pd

from tqdm import tqdm
from pprint import pformat
from Bio import SeqIO

from wisp_light.training.create_prediction import prediction
from wisp_light.training.utils import setup_logger, extract_majority_classification
from wisp_light.dataset.refSeqDataset import TAXO_LEVELS



class TaxoPredictor:
    def __init__(self, model_dir: str, params_file: str, predict_dir = None):
        self.model_dir = model_dir
        self.predict_dir = predict_dir
        self.temp_dir = './temp'

        if self.predict_dir:
            log_file = os.path.join( self.predict_dir, 'prediction.logs')
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
        else:
            log_file = None
        self.logger = setup_logger('TaxoPredictor', level=logging.INFO, log_file=log_file)

        # Load phylogenetic tree
        phylo_path = os.path.join(model_dir, "phylo_tree.txt")
        with open(phylo_path, 'rb') as jtree:
            self.phylo_tree = pickle.load(jtree)

        with open(params_file, 'r') as file:
            self.params = yaml.safe_load(file)

        self.model_dir = os.path.join(model_dir, "model")
        self.logger.info('Params for prediction\n' + pformat(self.params))
        self.logger.info(f"TaxoPredictor init for XGBoost models from directory: {model_dir}")
        self.logger.info("=" * 60)

        self.results_table = pd.DataFrame(columns=['fna_file', 'id_seq'] + TAXO_LEVELS )


    def predict_one_dir(self, fna_dir: str, raw_pred=False, verbose=False):

        if not os.path.isdir(fna_dir):
            raise NotADirectoryError(f"{fna_dir} n'est pas un répertoire valide.")
        fna_files = [f for f in os.listdir(fna_dir) if f.endswith(".fna")]
        self.logger.info(f"Prediction for all dir {fna_dir} contains {len(fna_files)} fna files.")

        if not fna_files:
            raise FileNotFoundError(f"Aucun fichier .fna trouvé dans {fna_dir}.")

        all_dir_results = {}

        for fna_file in tqdm(fna_files):
            fna_path = os.path.join(fna_dir, fna_file)
            if verbose:
                self.logger.info(f"Processing file: {fna_file}")
            fna_results = self.predict_one_fna(fna_path, raw_pred=raw_pred, verbose=verbose)
            all_dir_results.update(fna_results)

        return all_dir_results

    def predict_one_fna(self, fna_path: str, raw_pred=False, verbose=False):
        """Classify all sequences in a multi-fasta file (.fna)."""
        if not fna_path.endswith('.fna'):
            raise ValueError(f"Le fichier {fna_path} n'a pas une extension '.fna'")

        with open(fna_path, 'r', encoding='utf-8') as freader:
            genome_data = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')}

        all_results = {}

        if verbose:
            self.logger.info(f"Prediction for {os.path.basename(fna_path)} contains {len(genome_data)} sequences.")

        for id_sequence, dna_sequence in genome_data.items():
            try:
                raw_prediction = prediction(id_sequence, dna_sequence, self.params, self.phylo_tree, self.model_dir, self.temp_dir, self.logger)
                pred_taxons = extract_majority_classification(raw_prediction)

                if verbose:
                    self.logger.info(f"Raw prediction for genome id {id_sequence} of {os.path.basename(fna_path)}")
                    self.logger.info('\n' + pformat(raw_prediction))
                    self.logger.info(f"Taxo prediction for genome id {id_sequence} of {os.path.basename(fna_path)}")
                    self.logger.info('\n' + pformat( dict(pred_taxons)))

                raw_pred_to_level= {level: list(level_pred_dict.values()) for level, level_pred_dict in zip(TAXO_LEVELS, raw_prediction)}
                if raw_pred:
                    result_to_table = raw_pred_to_level
                    all_results[f"{os.path.basename(fna_path)}_id_{id_sequence}"] = raw_prediction

                else:
                    result_to_table =  pred_taxons
                    all_results[f"{os.path.basename(fna_path)}_id_{id_sequence}"] = pred_taxons


            except ValueError as e:
                self.logger.debug(f"Failed to predict for {os.path.basename(fna_path)} {id_sequence} : {str(e)}")
                result_to_table = {level: "Seq too short" for level in TAXO_LEVELS}
                all_results[f"{os.path.basename(fna_path)}_id_{id_sequence}"] = {"Seq too short": str(e)}

            row_data = {
                'fna_file': os.path.basename(fna_path),
                'id_seq': id_sequence,
                **{level: result_to_table.get(level, None) for level in TAXO_LEVELS}
            }
            self.results_table.loc[len(self.results_table)] = row_data

        return all_results


if __name__=='__main__':
    model_dir = "/home/hcourtei/Projects/MicroTaxo/codes/exp/model_base_complete_05_15_16_37"
    predict_dir = None #  "/home/hcourtei/Projects/MicroTaxo/codes/predict"
    params_file = "/home/hcourtei/Projects/MicroTaxo/codes/wisp_light/prediction/predict_params.yaml"
    fna_path = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data/GCF_000725405.1_ASM72540v1_genomic.fna"
    fna_dir = "/home/hcourtei/Projects/MicroTaxo/codes/predict/to_predict"

    # predictor = TaxoPredictor(model_dir, params_file, predict_dir=predict_dir)
    # all_results = predictor.predict_one_fna(fna_path, raw_pred=False)
    predictor = TaxoPredictor(model_dir, params_file, predict_dir=predict_dir)
    all_results = predictor.predict_one_fna(fna_path, raw_pred=True, verbose=True)
    # print(predictor.results_table.to_markdown())

    predictor = TaxoPredictor(model_dir, params_file, predict_dir=None)
    predictor.predict_one_dir(fna_dir, raw_pred=False, verbose=False)
    print(predictor.results_table.to_markdown())
