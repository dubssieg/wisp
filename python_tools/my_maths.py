from math import log


def entropy(explored_classes: int, probabilities: list[float]) -> float:
    """Computes the entropy for one given node

    Args:
        explored_classes (int): number of classes selected during exploration phase
        probabilities (list[float]): list of probabilities associated with explored classes

    Returns:
        float: entropy for node
    """
    if explored_classes > 1:
        return -sum([p*log(p, explored_classes) for p in probabilities])
    else:
        print("Not enough classes")
        return 0
