"Builds predictions from reads"
from os import remove
from collections import Counter
from xgboost import Booster, DMatrix
from numpy import argmax, amax, mean, ndarray
from treelib import Tree
from workspace.create_sample import build_sample


def make_prediction(
        model_path: str,
        parameters_path: str,
        datas_path: str,
        normalisation_func: str,
        read_identity_threshold: float
) -> list:
    """Does a prediction with a pre-calculated model

    Args:
        model_path (str): full path to model
        parameters_path (str): _description_
        datas_path (str): _description_
        normalisation_func (str): _description_
        read_identity_threshold (float): _description_

    Returns:
        list: _description_
    """

    # Creating booster object
    bst: Booster = Booster()
    bst.load_model(model_path)
    with open(parameters_path, "r", encoding='utf-8') as reader:
        bst.load_config(''.join([line for line in reader]))

    # Getting predictions
    predictions: ndarray = bst.predict(DMatrix(datas_path))

    return softmax(predictions, normalisation_func, read_identity_threshold)


def softmax(predictions: ndarray, func: str, reads_threshold: float) -> list:
    """Given a set of predictions, computes the consensus within it by ignoring some low-signifiance scores.

    Args:
        predictions (ndarray): array containing returns from the booster
        func (str): the function used to discriminate reads
        reads_threshold (float): between 0 and 1

    Raises:
        ValueError: if no read prediction is given (empty predictions array)

    Returns:
        list: a class for each read
    """

    if reads_threshold <= 0:
        return [argmax(a) for a in predictions]

    if reads_threshold > 1:
        reads_threshold = 1.0

    try:
        match func:
            case 'delta_mean':
                ret: list = [argmax(a) if amax(a)-mean(a) >
                             reads_threshold else False for a in predictions]
            case 'min_max':
                ret: list = [argmax(a) if min([amax(a)-p for p in a if p != amax(a)]) >
                             reads_threshold else False for a in predictions]
            case 'delta_sum':
                ret: list = [argmax(a) if amax(a) > (sum(a)-amax(a)) + reads_threshold
                             else False for a in predictions]
            case _:
                ret: list = [argmax(a) for a in predictions]
    except ValueError:
        ret: list = [argmax(a) for a in predictions]
    if len(ret) == 0:
        raise ValueError(
            "There's no read to evaluate, your entry data might be broken.")
    elif not all(not p for p in ret):
        # print(f"All reads with a threshold inferior at {round(reads_threshold,ndigits=2)} for function {func} have been purged.")
        return ret
    else:
        return softmax(predictions, func, reads_threshold-0.05)


def prediction(id_sequence: str, dna_sequence: str, params: dict, tree: Tree, threshold: float, parameters: str) -> str | list:
    """Creates a prediction for a read.

    Args:
        id_sequence (str): identifier for sequence (fasta header)
        dna_sequence (str): the full DNA sequence
        params (dict): params global dict
        tree (Tree): taxonomy tree built before
        threshold (float): threshold to consider a taxa as accurate
        parameters (str): path to parameters file

    Returns:
        str | list: prediction results
    """
    if len(dna_sequence) >= params['read_size']:
        # Creating the sample
        sample_output_path: str = build_sample(
            parameters,
            dna_sequence,
            id_sequence
        )
    else:
        return 'REJECTED'

    # Evaluate at one level
    results: list[dict] = [{} for _ in range(5)]
    kept_taxas: list = ['Root']
    for i, _ in enumerate(['root', 'domain', 'phylum', 'group', 'order'], start=0):
        mappings_taxa: dict = {
            node.data.code: node.tag for node in tree.filter_nodes(lambda x: tree.depth(x) == i+1)}
        # Use the tree to select next level
        for taxa in tree.filter_nodes(lambda x: tree.depth(x) == i):
            if taxa.tag in kept_taxas:
                # print(f"Processing {taxa.tag} @ {level}")
                results[i][taxa.tag] = {mappings_taxa[key]: value for key, value in Counter(make_prediction(
                    taxa.data.model_path,
                    taxa.data.config_path,
                    sample_output_path,
                    normalisation_func='delta_mean',
                    read_identity_threshold=0.8
                )).items() if key is not False}
        kept_taxas = [lower_taxa for counter in results[i].values(
        ) for lower_taxa, count in counter.items() if count > threshold*sum(list(counter.values()))]

    # Cleaning temp files
    remove(sample_output_path)

    return results
