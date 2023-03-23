from xgboost import Booster, DMatrix
from numpy import argmax, amax, mean, ndarray


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
    "Given a set of predictions, computes the consensus within it by ignoring some low-signifiance scores."

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
        print(
            f"All reads with a threshold inferior at {reads_threshold} for function {func} have been purged.")
        return ret
    else:
        return softmax(predictions, func, reads_threshold-0.05)
