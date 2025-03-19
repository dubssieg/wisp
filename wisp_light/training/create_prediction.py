import time
import os
import sys
from collections import Counter
from utils import softmax
from xgboost import Booster, DMatrix
from xgboost.core import XGBoostError
from create_database import encoder,splitting, counter_kmer

sys.path.append('../../..')
from wisp.wisp_light.dataset.refSeqDataset import TAXO_LEVELS

# LEVELS = ['root', 'domain', 'phylum', 'group', 'order']

def make_prediction(model_path, datas_path, normalisation_func,read_identity_threshold) -> list:
    """Does a prediction with a pre-calculated model
    """
    # Creating booster object
    bst = Booster()
    bst.load_model(model_path)
    config_path = model_path.replace('.json', '_params.json')

    with open(config_path, "r", encoding='utf-8') as reader:
        bst.load_config(''.join([line for line in reader]))
    try :
        predictions = bst.predict(DMatrix(datas_path+"?format=libsvm")) #ndarray
        # print(predictions)
    except XGBoostError as e:
        raise RuntimeError(f"Error processing {datas_path}: {e}") from e  # Propager l'erreur avec plus de contexte        return None

    list_resul = softmax(predictions, normalisation_func, read_identity_threshold)
    # print(list_resul)
    return list_resul



def build_sample(params: dict, dna_sequence: str, id_sequence: str, sample_output_path) -> None :
    "Builds a json file with taxa levels as dict information"
    os.makedirs(os.path.dirname(sample_output_path), exist_ok=True)     # Writing the database

    my_encoder: dict = encoder(ksize=params['ksize'])
    all_reads = splitting(dna_sequence.upper(), params['read_size'], params['max_sampling'], shift_ratio=params['shift_ratio'])

    with open(sample_output_path, 'w', encoding='utf-8') as jdb:
        # Counting kmers inside each read
        counters =  [counter_kmer(read, params['pattern']) for read in all_reads]
        encoded: list = [{my_encoder[k]:v for k, v in cts.items()} for cts in counters]   # Encoding reads for XGBoost

        for sample in encoded:
            jdb.write(f"0 {' '.join([str(k)+':'+str(v) for k,v in sample.items()])} #{id_sequence}\n")
            # Each read is a dict with code:count for kmer

def prediction(id_sequence: str, dna_sequence: str, params: dict, tree, model_dir:str, val_dir) ->list:
    """Creates a prediction for a read.

    Args:
        id_sequence (str): identifier for sequence (fasta header)
        dna_sequence (str): the full DNA sequence
        params (dict): params global dict
        tree (Tree): taxonomy tree built before

    Returns:
        list: prediction results
    """
    if len(dna_sequence) < params['read_size']:
        raise ValueError(f"DNA sequence too short {len(dna_sequence):,} and minimum required: {params['read_size']:,}")

    file = f"unk_sample_{str(time.time()).replace('.', '_')}_{id_sequence.replace(' ', '_')}.txt"
    sample_output_path = f"{val_dir}/temp/{file}"

    build_sample(params, dna_sequence, id_sequence, sample_output_path)
    # Evaluate at one level
    results: list[dict] = [{} for _ in range(5)]
    kept_taxas = ['Root']
    for id_level, level in enumerate(['root'] + TAXO_LEVELS[:-1]):
        mappings_taxa = {node.data.code: node.tag for node in tree.filter_nodes(lambda x: tree.depth(x) == id_level+1)}
        # Use the tree to select next level
        for taxa in tree.filter_nodes(lambda x: tree.depth(x) == id_level):
            if taxa.tag not in kept_taxas:
                continue  # On ignore les taxons non sélectionnés
            if taxa.data.model_path is None:
                print(f"⚠️ Avertissement: Aucun modèle pour {taxa.tag} (niveau {level})")
                continue

            model_path = os.path.join(model_dir, taxa.data.model_path)
            predictions = make_prediction(model_path, sample_output_path,
                                          normalisation_func=params['normalisation_func'],
                                          read_identity_threshold=params['read_identity_threshold'])
            # print(f"level{i} predictions len {len(predictions)} \n {predictions}")
            ## TODO enlever le passage en liste dans softmax et traiter différement les False
            count_items = Counter(predictions).items()
            results[id_level][taxa.tag] = {mappings_taxa[key]: value for key, value in count_items if key is not False}

        # Sélection des taxons à garder pour le prochain niveau
        kept_taxas = [lower_taxa for counter in results[id_level].values() for lower_taxa, count in counter.items()
            if count > params['threshold'] * sum(counter.values())]

    os.remove(sample_output_path)
    return results

