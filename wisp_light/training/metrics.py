from collections import defaultdict
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd
import pickle


from wisp_light.dataset.refSeqDataset import TAXO_LEVELS


class ConfusionMatrixTracker:
    """
    Accumule et génère des matrices de confusion pour chaque niveau taxonomique.

    Attributes
    ----------
    true_labels : dict of list
       Listes des labels vrais par niveau.
    pred_labels : dict of list
       Listes des labels prédits par niveau.
    unknown_pred : str
       Valeur utilisée pour les prédictions manquantes.
    taxonomy_df : pandas.DataFrame or None
       DataFrame des hiérarchies taxonomiques uniques (après build_taxonomy_df).
    """

    def __init__(self):
        """
        Initialise les stockages pour labels vrais et prédits.

        Définit des listes vides pour chaque niveau de TAXO_LEVELS et
        configure la valeur par défaut pour les prédictions inconnues.
        """
        # Stocke les vraies et prédictions pour chaque niveau
        self.true_labels = defaultdict(list)
        self.pred_labels = defaultdict(list)
        self.unknown_pred = "unknown"


    def update(self, true_labels, pred_labels):
        """
        Met à jour les listes de labels vrais et prédits pour chaque niveau.

        Parameters
        ----------
        true_labels : dict
            Dictionnaire {niveau: label_vrai}. Aucun niveau ne doit être None.
        pred_labels : dict
            Dictionnaire {niveau: label_prédit}. Les niveaux absents ou None
            seront remplacés par `self.unknown_pred`.
        """
        for level in TAXO_LEVELS:
            true_value = true_labels.get(level)
            assert true_value is not None, f"True label for level '{level}' is None. All labels {true_labels}"
            pred_value = pred_labels.get(level, self.unknown_pred)  # Remplace None par "Unknown"

            if pred_value is None:
                pred_value = self.unknown_pred  # Remplace explicitement les None

            self.true_labels[level].append(true_value)
            self.pred_labels[level].append(pred_value)

    def build_taxonomy_df(self):
        """
        Construit et stocke un DataFrame des taxonomies uniques.

        Assemble les labels vrais accumulés en un array NumPy de forme
        (n_samples, n_levels), puis crée un DataFrame trié et sans doublons.

        After calling this method, `self.taxonomy_df` is a DataFrame with one
        row per unique taxonomy path, sorted lexicographically by TAXO_LEVELS.
        The tracker is also pickled to "conf_matrix_tracker.pkl".
        """
        data = []

        # Itérer sur les indices des échantillons
        num_samples = len(self.true_labels[TAXO_LEVELS[0]])  # Nombre d'échantillons
        for i in range(num_samples):
            row = {level: self.true_labels[level][i] for level in TAXO_LEVELS}  # Utiliser l'indice 'i'
            data.append(row)
        taxonomy_df = pd.DataFrame(data)
        taxonomy_df.sort_values(by=TAXO_LEVELS, inplace=True)
        taxonomy_df = taxonomy_df.drop_duplicates().reset_index(drop=True)
        # print("ehehe", taxonomy_df.to_markdown())
        # print("after sorting \n", taxonomy_df)
        self.taxonomy_df = taxonomy_df
        with open("conf_matrix_tracker.pkl", "wb") as f:
            pickle.dump(self, f)

    def get_confusion_matrix(self, level):
        """
        Calcule la matrice de confusion pour un niveau taxonomique donné.

        Parameters
        ----------
        level : str
            Nom d’un niveau dans TAXO_LEVELS (ex. 'phylum', 'class', …).

        Returns
        -------
        pd.DataFrame
            Matrice de confusion (labels vrais en lignes, prédits en colonnes).
            Les classes sont ordonnées selon `self.taxonomy_df[level]`, avec
            `self.unknown_pred` ajoutée si nécessaire.

        Raises
        ------
        AttributeError
            Si `build_taxonomy_df` n’a pas encore été appelé.
        """

        if not hasattr(self, 'taxonomy_df'):
            raise AttributeError("Attributes 'taxonomy_df' not present , call before self.build_taxonomy_df.")

        all_y_true = self.true_labels[level]
        all_y_pred = self.pred_labels[level]

        sorted_classes = self.taxonomy_df[level].drop_duplicates().tolist()

        if self.unknown_pred in all_y_pred:
            sorted_classes.append(self.unknown_pred)

        if not sorted_classes:
            return pd.DataFrame()

        matrix = confusion_matrix(all_y_true, all_y_pred, labels=sorted_classes)
        confusion_df = pd.DataFrame(matrix, index=sorted_classes, columns=sorted_classes)

        return confusion_df


    def calculate_separator_indices(self, level_marker='phylum', level_index='family'):
        """
        Identifie les indices où le label de `level_marker` change en triant
        par `level_index`.

        Parameters
        ----------
        level_marker : str
            Niveau supérieur dont on suit les changements (ex. 'phylum').
        level_index : str
            Niveau inférieur servant d’index pour le tri (ex. 'family').

        Returns
        -------
        list of int
            Positions dans la liste triée où la valeur de `level_marker` change.

        Raises
        ------
        AttributeError
            Si `build_taxonomy_df` n’a pas encore été appelé.
        """
        if not hasattr(self, 'taxonomy_df'):
            raise AttributeError("L'attribut 'taxonomy_df' not present , call before self.build_taxonomy_df")
        # unique_df = self.taxonomy_df.drop_duplicates(subset=[level_2])
        #
        # # Trier le DataFrame selon le niveau inférieur

        # sorted_df = unique_df.sort_values(by=[level_2]).reset_index(drop=True)
        unique_df = self.taxonomy_df.drop_duplicates(subset=[level_index]).reset_index()
        separator_indices = []

        previous_level_1_value = None
        for i, row in unique_df.iterrows():
            current_level_1_value = row[level_marker]
            if previous_level_1_value is not None and current_level_1_value != previous_level_1_value:
                separator_indices.append(i)  # Ajouter l'indice où level_1 change

            previous_level_1_value = current_level_1_value

        return separator_indices


def compute_accuracy_from_conf_matrix_df(conf_mat_level):
    """
    Calcule l’exactitude globale à partir d’une matrice de confusion.

    Parameters
    ----------
    conf_mat_level : pd.DataFrame
        Matrice de confusion (carrée) pour un niveau donné.

    Returns
    -------
    float
        Taux de bonnes classifications = trace / somme de tous les éléments.
        Retourne 0.0 si la somme totale est nulle.
    """
    correct_predictions = np.diag(conf_mat_level).sum()
    total_predictions = conf_mat_level.sum().sum()
    accuracy = correct_predictions / total_predictions if total_predictions > 0 else 0.0
    return accuracy

if __name__ == "__main__":
    from wisp_light.visu.plots_tools import plot_conf_mat

    with open("conf_matrix_tracker.pkl", "rb") as f:
        tracker = pickle.load(f)

    # Tu peux maintenant utiliser toutes les méthodes et attributs :
    print(tracker.taxonomy_df.to_markdown())  # Voir le DataFrame taxonomique
    for level in  ['phylum', 'class', 'order', 'family']:
        level_marker = 'phylum'
        conf_mat_level = tracker.get_confusion_matrix(level)
        print(conf_mat_level)  # Voir la matrice de confusion du niveau "phylum"
        # print(tracker.true_labels["phylum"][:5])  # Voir quelques labels
        # print(tracker.pred_labels["phylum"][:5])  # Voir quelques prédictions
        if level != level_marker:
            separator_indices = tracker.calculate_separator_indices(level_marker=level_marker, level_index=level)
        else:
            separator_indices = None
        print(separator_indices)
        plot_conf_mat(conf_mat_level, level, separator_indices, filename=None)