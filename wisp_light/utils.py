import json
import logging

def setup_logger(name: str, level=logging.INFO, console_level=logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.StreamHandler()
    # handler.setLevel(console_level)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%H:%M')) #''%Y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
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