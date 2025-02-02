import json
import logging
from numpy import argmax, amax, mean, ndarray

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

    reads_threshold = min(1.0, reads_threshold)  # Limite supérieure à 1.0

    try:
        match func:
            case 'delta_mean':
                ret = [argmax(a) if amax(a)-mean(a) >
                             reads_threshold else False for a in predictions]
            case 'min_max':
                ret = [argmax(a) if min([amax(a)-p for p in a if p != amax(a)]) >
                             reads_threshold else False for a in predictions]
            case 'delta_sum':
                ret = [argmax(a) if amax(a) > (sum(a)-amax(a)) + reads_threshold
                             else False for a in predictions]
            case _:
                ret = [argmax(a) for a in predictions]
    except ValueError:
        ret = [argmax(a) for a in predictions]
    if len(ret) == 0:
        raise ValueError("There's no read to evaluate, your entry data might be broken.")
    elif not all(not p for p in ret):
        # print(f"All reads with a threshold inferior at {round(reads_threshold,ndigits=2)} for function {func} have been purged.")
        return ret
    else:
        return softmax(predictions, func, reads_threshold-0.05)


def setup_logger(name: str, level=logging.INFO, log_file=None) -> logging.Logger:
    format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%H:%M')
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.StreamHandler()
    # handler.setLevel(console_level)
    handler.setFormatter(format)#''%Y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(format)  # ''%Y-%m-%d %H:%M:%S'))
        logger.addHandler(file_handler)

    return logger

def display_json_preview(file_path: str, num_elements: int = 5):
    """Affiche un aperçu formaté d'un fichier JSON.

    Args:
        file_path (str): Chemin du fichier JSON.
        num_elements (int): Nombre d'éléments à afficher.
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            if isinstance(data, dict):  # Si c'est un objet JSON
                preview = {key: data[key] for key in list(data.keys())[:num_elements]}
            elif isinstance(data, list):  # Si c'est une liste JSON
                preview = data[:num_elements]
            else:
                print("Le fichier JSON contient un format inattendu.")
                return

            print(json.dumps(preview, indent=4))  # Affichage formaté
    except FileNotFoundError:
        print(f"Le fichier {file_path} n'existe pas.")
    except json.JSONDecodeError:
        print(f"Erreur de décodage JSON dans le fichier {file_path}.")
    except Exception as e:
        print(f"Une erreur s'est produite : {e}")