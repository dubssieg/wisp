import json
import logging
import os
import psutil
import sys

from numpy import argmax, amax, mean, ndarray,array, vectorize
from wisp_light.dataset.refSeqDataset import TAXO_LEVELS


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
                ret = [argmax(a) if amax(a)-mean(a) > reads_threshold else False for a in predictions]
            case 'min_max':
                ret = [argmax(a) if min([amax(a)-p for p in a if p != amax(a)]) > reads_threshold else False for a in predictions]
            case 'delta_sum':
                ret = [argmax(a) if amax(a) > (sum(a)-amax(a)) + reads_threshold else False for a in predictions]
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
        # print(f"one more softmax with reads_threshold {reads_threshold-0.05}")
        return softmax(predictions, func, reads_threshold-0.05)

#
# def softmax2(predictions: ndarray, func: str, reads_threshold: float,min_threshold: float = 0) -> ndarray:
#     """Given a set of predictions, computes the consensus within it by ignoring some low-signifiance scores."""
#
#     if reads_threshold <= min_threshold:
#         return predictions  # Retourne directement le tableau d'origine si le seuil est <= 0
#
#     reads_threshold = min(1.0, reads_threshold)  # Limite supérieure à 1.0
#
#     def evaluate_read(a):
#         """Applique la logique de discrimination selon le type de fonction."""
#         if func == 'delta_mean' and (amax(a) - mean(a)) > reads_threshold:
#             return a
#         elif func == 'min_max' and (amax(a) - min([p for p in a if p != amax(a)])) > reads_threshold:
#             return a
#         elif func == 'delta_sum' and amax(a) > (sum(a) - amax(a)) + reads_threshold:
#             return a
#         return None  # Retourne None pour les prédictions à ignorer.
#
#     # Applique la logique de discrimination sur chaque prédiction
#     vectorized_evaluation = vectorize(evaluate_read, otypes=[ndarray])  # Otypes assure que le retour est un array numpy
#     ret = vectorized_evaluation(predictions)
#
#     # Filtre les valeurs None et retourne un tableau numpy
#     ret = ret[ret is not None]
#
#     # Vérifie si des prédictions valides existent
#     if ret.size == 0:
#         if reads_threshold > min_threshold:
#             # Si aucune prédiction n'est valide, on diminue le seuil et réessaye
#             return softmax2(predictions, func, reads_threshold - 0.05, min_threshold)
#         else:
#             # Si le seuil a atteint la limite minimale et aucune prédiction n'est valide, on lève une exception
#             raise ValueError("No valid predictions, data might be broken.")
#
#
#     return ret



def log_resource_usage(logger: logging.Logger, label: str = ""):
    pid = os.getpid()
    process = psutil.Process(pid)
    num_threads = process.num_threads()
    num_children = len(process.children(recursive=True))
    memory = process.memory_info().rss / 1024 / 1024  # in MB
    cpu_percent = process.cpu_percent(interval=0.1)

    green = "\033[92m"
    reset = "\033[0m"
    msg = (
        f"{green}[{label}] PID={pid} 🧵 Threads: {num_threads} | 👶 Subprocesses: {num_children} | "
        f"🧠 RAM: {memory:.2f} MB | 🔥 CPU%: {cpu_percent}{reset}"
    )

    logger.info(msg)


def extract_majority_classification(sequence):
    """
    Extrait la classification majoritaire pour chaque niveau (domain, phylum, group, order, family) d'une séquence donnée.
    Args:
        sequence (list): Liste de dictionnaires représentant les niveaux taxonomiques avec leurs pourcentages.
    Returns:
        dict: Dictionnaire avec la classification majoritaire pour chaque niveau.
    """
    classification = dict(zip(TAXO_LEVELS, [None] * len(TAXO_LEVELS)))

    # Niveau de classification pour chaque taxon (ordre défini)

    # Pour chaque niveau, on prend le taxon majoritaire basé sur le pourcentage
    for id_level, level in enumerate(sequence):
        if not level:  # Si le niveau est vide, on passe
            continue
        for key, value in level.items():
            if isinstance(value, dict):  # Si la valeur est un sous-dictionnaire
                major_taxon = max(value.items(), key=lambda x: x[1])[0]
                # Assigner le taxon majoritaire selon l'ordre des clés
                if classification[TAXO_LEVELS[id_level]] is None:
                    classification[TAXO_LEVELS[id_level]] = major_taxon

    return classification



def setup_logger(name: str, level=logging.INFO, log_file=None) -> logging.Logger:
    format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%H:%M')
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stdout)
    # handler.setLevel(console_level)
    handler.setFormatter(format)#''%Y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)

    if log_file:
        sys.stdout.write(log_file + "\n")
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