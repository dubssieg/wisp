import time
import os
from collections import Counter
from xgboost import Booster, DMatrix
from xgboost.core import XGBoostError

from .utils import softmax
from .create_database import encoder,splitting, counter_kmer
from wisp_light.dataset.refSeqDataset import TAXO_LEVELS


def make_prediction(model_path, datas_path, normalisation_func,read_identity_threshold) -> list:
    """
    Makes predictions using a pre-trained XGBoost model.

    Parameters
    ----------
    model_path : str
        Path to the saved XGBoost model file.
    datas_path : str
        Path to the input dataset in LibSVM format.
    normalisation_func : function
        The normalisation function applied to the predictions.
    read_identity_threshold : float
        Threshold for read identity.

    Returns
    -------
    list
        List of predictions after applying the softmax and normalisation function.
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
        with open(datas_path, "r", encoding="utf-8") as f:
            for idx, line in enumerate(f):
                # on splitpe la ligne en tokens, on ignore le label (premier token)
                parts = line.strip().split()
                # on collecte tous les indices de features
                features = {int(tok.split(":")[0]) for tok in parts[1:] if ":" in tok}
                # si 0 est parmi les clés…
                if 0 in features:
                    print(f"Ligne {idx} contient la feature 0 → {line.strip()}")
                if max(features) >= bst.num_features():
                    print(f"Ligne {idx} contient un index trop grand (max={max(features)}) → {line.strip()}")

        raise RuntimeError(f"Error processing {datas_path}: {e}") from e  # Propager l'erreur avec plus de contexte        return None

    list_resul = softmax(predictions, normalisation_func, read_identity_threshold)
    # print(list_resul)
    return list_resul



def build_sample(params: dict, dna_sequence: str, id_sequence: str, sample_output_path) -> None :
    """
    Builds a JSON file representing k-mer counts for a given DNA sequence.

    Parameters
    ----------
    params : dict
        Dictionary of parameters, including k-mer size and read size.
    dna_sequence : str
        The full DNA sequence.
    id_sequence : str
        The identifier for the sequence (e.g., a FASTA header).
    sample_output_path : str
        The path where the resulting sample file will be saved.

    Returns
    -------
    None
        This function does not return any value.
    """
    os.makedirs(os.path.dirname(sample_output_path), exist_ok=True)     # Writing the database

    my_encoder: dict = encoder(ksize=params['ksize'])
    # num_features = len(set(my_encoder.values()))
    # print(f"Nombre total de features dans build_sample : {num_features}")
    all_reads = splitting(dna_sequence.upper(), params['read_size'], params['max_sampling'], shift_ratio=params['shift_ratio'])

    with open(sample_output_path, 'w', encoding='utf-8') as jdb:
        # Counting kmers inside each read
        counters =  [counter_kmer(read, params['pattern']) for read in all_reads]
        encoded: list = [{my_encoder[k]:v for k, v in cts.items()} for cts in counters]   # Encoding reads for XGBoosts]
        max_encode = int('3'*params['ksize'])

        for sample in encoded:
            sample.pop(max_encode, None)

            feats = [f"{k}:{v}" for k, v in sample.items()]
            jdb.write(f"0 {' '.join(feats)} #{id_sequence}\n")


def prediction(id_sequence: str, dna_sequence: str, params: dict, tree, model_dir, val_dir, logger) ->list:
    """
    Creates a prediction for a given DNA sequence using a taxonomy tree and pre-trained models.

    Parameters
    ----------
    id_sequence : str
        The identifier for the sequence (e.g., a FASTA header).
    dna_sequence : str
        The full DNA sequence.
    params : dict
        Dictionary of parameters, including read size and threshold values.
    tree : Tree
        The taxonomy tree used for hierarchical predictions.
    model_dir : str
        Directory containing pre-trained models for each taxonomic level.
    val_dir : str
        Directory for storing temporary files during prediction.
    logger : Logger
        Logger object for logging information and warnings.

    Returns
    -------
    list
        A list of dictionaries containing the predictions for each taxonomic level.
    """
    if len(dna_sequence) < params['read_size']:
        raise ValueError(f"DNA sequence too short {len(dna_sequence):,} and minimum required: {params['read_size']:,}")
    # if len(dna_sequence) < params['read_size']:
        # Catch the exception and return None (skip this sequence)
        # logger.warning(f"⚠️ Sequence '{id_sequence}' is too short ({len(dna_sequence)}) and will be skipped.")
        # return None  # Skip this sequence
    file = f"unk_sample_{str(time.time()).replace('.', '_')}_{id_sequence.replace(' ', '_')}.txt"
    sample_output_path = f"{val_dir}/temp/{file}"
    build_sample(params, dna_sequence, id_sequence, sample_output_path)
    # Evaluate at one level
    results: list[dict] = [{} for _ in range(len(TAXO_LEVELS))]
    kept_taxas = ['Root']

    for id_level, level in enumerate(['root'] + TAXO_LEVELS[:-1]):
        mappings_taxa = {node.data.code: node.tag for node in tree.filter_nodes(lambda x: tree.depth(x) == id_level+1)}
        # Use the tree to select next level
        for taxa in tree.filter_nodes(lambda x: tree.depth(x) == id_level):
            if taxa.tag not in kept_taxas:
                continue  # On ignore les taxons non sélectionnés
            if taxa.data.model_path is None:
                logger.info(f"⚠️ Avertissement: Aucun modèle pour {taxa.tag} (niveau {level}) taxa.tag {taxa.tag}")
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

