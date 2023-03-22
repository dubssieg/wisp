from xgboost import Booster, DMatrix
from numpy import argmax, amax, mean, ndarray


def predict(model: Booster, datas_path: str, limit: float, func: str) -> list:

    # Getting predictions
    predictions: ndarray = model.predict(DMatrix(datas_path))

    return softmax(predictions, func, limit)


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
