import os
import pickle
import logging
import yaml
import pandas as pd

from tqdm import tqdm
from pprint import pformat
from Bio import SeqIO
from datetime import datetime

import sys
sys.path.append('../..')
from wisp_light.training.create_prediction import prediction
from wisp_light.training.utils import setup_logger, extract_majority_classification
from wisp_light.dataset.refSeqDataset import TAXO_LEVELS



class TaxoPredictor:
    """
    Classe permettant de prédire la classification taxonomique à partir de fichiers génomiques `.fna`
    à l'aide d'un modèle XGBoost entraîné.

    Cette classe charge un arbre phylogénétique, des paramètres de prédiction depuis un fichier YAML,
    et utilise un modèle sauvegardé pour effectuer des prédictions sur un ou plusieurs fichiers `.fna`.

    Parameters
    ----------
    model_dir : str
        Chemin vers le répertoire contenant le modèle entraîné et l’arbre phylogénétique (`phylo_tree.txt`).
    params_file : str
        Chemin vers le fichier YAML contenant les paramètres nécessaires à la prédiction.
    predict_dir : str, optional
        Répertoire où seront sauvegardés les fichiers de log et les résultats. Si None, les logs ne seront pas enregistrés sur disque.

    Attributes
    ----------
    model_dir : str
        Chemin vers le sous-dossier contenant les fichiers du modèle XGBoost.
    predict_dir : str or None
        Répertoire où les résultats et les logs sont sauvegardés.
    temp_dir : str
        Répertoire temporaire utilisé lors des étapes de prédiction.
    logger : logging.Logger
        Logger utilisé pour enregistrer les messages d'information et de débogage.
    phylo_tree : dict
        Arbre phylogénétique chargé depuis un fichier `pickle`, utilisé pour guider les prédictions.
    results_table : pandas.DataFrame
        Table contenant les résultats de prédiction pour chaque séquence analysée.
    """
    def __init__(self, model_dir: str, params_file: str, predict_dir:str):
        self.model_dir = model_dir
        self.predict_dir = predict_dir
        self.temp_dir = './temp'

        if self.predict_dir:
            os.makedirs(self.predict_dir, exist_ok=True)
            base_log_path = os.path.join(self.predict_dir, 'prediction.logs')
            log_file = base_log_path
            index = 1
            while os.path.exists(log_file):
                log_file = os.path.join(self.predict_dir, f"prediction_{index}.logs")
                index += 1
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
        """
        Prédit la taxonomie pour tous les fichiers `.fna` présents dans un répertoire.

        Parameters
        ----------
        fna_dir : str
            Chemin vers un répertoire contenant des fichiers `.fna`.
        raw_pred : bool, default=False
            Si True, retourne les prédictions brutes pour chaque niveau taxonomique.
            Si False, retourne les prédictions finales basées sur un vote majoritaire.
        verbose : bool, default=False
            Si True, affiche des informations détaillées pendant le traitement.

        Returns
        -------
        dict
            Dictionnaire associant les identifiants de séquences aux résultats de prédiction.

        Raises
        ------
         NotADirectoryError
            Levée si `fna_dir` n'est pas un répertoire valide.
        FileNotFoundError
            Levée si aucun fichier `.fna` n'est trouvé dans le répertoire.
        """

        if not os.path.isdir(fna_dir):
            raise NotADirectoryError(f"{fna_dir} n'est pas un répertoire valide.")
        fna_files = [f for f in os.listdir(fna_dir) if f.endswith(".fna")]
        self.logger.info(f"Prediction for whole dir {fna_dir},  containing {len(fna_files)} fna files.")

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
        """
        Prédit la taxonomie pour chaque séquence contenue dans un fichier multi-FASTA (.fna).

        Parameters
        ----------
        fna_path : str
            Chemin vers un fichier `.fna` contenant une ou plusieurs séquences génomiques.
        raw_pred : bool, default=False
             Si True, retourne les prédictions brutes pour chaque niveau taxonomique.
            Si False, retourne les prédictions finales issues d'un vote majoritaire.
        verbose : bool, default=False
            Si True, affiche des informations détaillées pendant le traitement.

        Returns
        -------
        dict
            Dictionnaire associant les identifiants de séquence (avec référence au fichier) aux résultats de prédiction.

        Raises
        ------
        ValueError
            Levée si le fichier fourni ne possède pas l'extension `.fna`.
        """
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

    def save_result_to_csv(self):

        """
        Sauvegarde la table des résultats dans un fichier CSV avec un nom horodaté.

        Si un fichier avec le même nom existe déjà, un indice (_1, _2, ...) est ajouté
        au nom pour éviter l'écrasement.

        Le format de l’horodatage est `_JJ_MM_HH_MM`. Le fichier CSV est sauvegardé dans
        le répertoire spécifié par `self.predict_dir`. Si `self.predict_dir` n’est pas défini,
        un message est enregistré dans les logs et aucune sauvegarde n’est effectuée.

        Parameters
        ----------
            None

        Returns
        --------
            None
        """
        if not self.predict_dir:
            self.logger.info("predict_dir is not set — résultats non sauvegardés.")
            return

        os.makedirs(self.predict_dir, exist_ok=True)
        date_suffix = datetime.now().strftime("_%d_%m_%H_%M")
        base_name = f"results_table{date_suffix}"
        file_path = os.path.join(self.predict_dir, f"{base_name}.csv")

        index = 1
        while os.path.exists(file_path):
            file_path = os.path.join(self.predict_dir, f"{base_name}_{index}.csv")
            index += 1

        self.results_table.to_csv(file_path, index=False, sep=";")
        self.logger.info(f"Table des résultats sauvegardée dans : {file_path}")


if __name__=='__main__':
    model_dir = "/home/hcourtei/Projects/MicroTaxo/codes/exp/model_base_complete_05_15_16_37"
    predict_dir = "/home/hcourtei/Projects/MicroTaxo/codes/predict/predicted"  # permet de sauvegarder les log et les predictions
    params_file = "/home/hcourtei/Projects/MicroTaxo/codes/wisp_light/prediction/predict_params.yaml"
    fna_path = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data/GCF_000725405.1_ASM72540v1_genomic.fna"
    fna_dir = "/home/hcourtei/Projects/MicroTaxo/codes/predict/to_predict"



    predictor = TaxoPredictor(model_dir, params_file,predict_dir=predict_dir)
    all_results = predictor.predict_one_fna(fna_path, raw_pred=True, verbose=True)
    print(predictor.results_table.to_markdown())

    # predictor.save_result_to_csv()
    #
    # predictor = TaxoPredictor(model_dir, params_file, predict_dir=predict_dir)
    # predictor.predict_one_dir(fna_dir, raw_pred=False, verbose=False)
    # predictor.save_result_to_csv()

