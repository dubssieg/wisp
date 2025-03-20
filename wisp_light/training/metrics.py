from collections import defaultdict
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd
import sys
sys.path.append('..')
from wisp.wisp_light.dataset.refSeqDataset import TAXO_LEVELS


class ConfusionMatrixTracker:
    def __init__(self):
        # Stocke les vraies et prédictions pour chaque niveau
        self.true_labels = defaultdict(list)
        self.pred_labels = defaultdict(list)


    def update(self, true_labels, pred_labels):
        """
        Met à jour les listes de vraies étiquettes et de prédictions.
        Args:
            true_labels (dict): {niveau: classe_vraie}
            pred_labels (dict): {niveau: classe_prédite}
        """
        for level in TAXO_LEVELS:
            true_value = true_labels.get(level)
            pred_value = pred_labels.get(level)

            if true_value is not None and pred_value is not None:  # Exclure les None
                self.true_labels[level].append(true_value)
                self.pred_labels[level].append(pred_value)

    def build_taxonomy_df(self):
        """
        Crée un DataFrame avec toutes les hiérarchies taxonomiques pour chaque échantillon
        en utilisant uniquement les labels de vérité terrain (ground truth).
        """
        data = []

        # Itérer sur les indices des échantillons
        num_samples = len(self.true_labels[TAXO_LEVELS[0]])  # Nombre d'échantillons
        for i in range(num_samples):
            row = {level: self.true_labels[level][i] for level in TAXO_LEVELS}  # Utiliser l'indice 'i'
            data.append(row)
        taxonomy_df = pd.DataFrame(data)
        # print("before sorting \n", taxonomy_df)
        taxonomy_df.sort_values(by=TAXO_LEVELS, inplace=True)
        taxonomy_df = taxonomy_df.drop_duplicates().reset_index(drop=True)
        # print("after sorting \n", taxonomy_df)
        self.taxonomy_df = taxonomy_df

    def get_confusion_matrix(self, level):
        """
        Retourne une matrice de confusion pour un niveau donné.

        Args:
            level (str): Un des niveaux taxonomiques ('domain', 'phylum', etc.)

        Returns:
            pd.DataFrame: Matrice de confusion avec noms des classes.
        """
        if level not in self.true_labels:
            raise ValueError(f"Niveau invalide: {level}")

        if not hasattr(self, 'taxonomy_df'):
            raise AttributeError("L'attribut 'taxonomy_df' not present , call before self.build_taxonomy_df.")

        y_true = self.true_labels[level]
        y_pred = self.pred_labels[level]

        # Filtrer les None
        filtered_true = [label for label in y_true if label is not None]
        filtered_pred = [label for label in y_pred if label is not None]

        sorted_classes = self.taxonomy_df[level].drop_duplicates().tolist()
        if not sorted_classes:
            return pd.DataFrame()

        matrix = confusion_matrix(filtered_true, filtered_pred, labels=sorted_classes)

        return pd.DataFrame(matrix, index=sorted_classes, columns=sorted_classes)


    def calculate_separator_indices(self, level_1='phylum', level_2='family'):
        """
        Calcule les indices où un changement du niveau level_1 (ex: phylum) se produit
        en regardant l'ordre du niveau level_2 (ex: family).
        Args:
            level_1 (str): Niveau supérieur (ex: 'phylum')
            level_2 (str): Niveau inférieur servant d'indexation (ex: 'family')

        Returns:
            separator_indices (list): Liste des indices où level_1 change selon level_2
        """
        if not hasattr(self, 'taxonomy_df'):
            raise AttributeError("L'attribut 'taxonomy_df' not present , call before self.build_taxonomy_df")
        separator_indices = []
        # Trier le DataFrame selon le niveau inférieur
        sorted_df = self.taxonomy_df.sort_values(by=[level_2]).reset_index(drop=True)
        previous_level_1_value = None
        for i, row in sorted_df.iterrows():
            current_level_1_value = row[level_1]
            if previous_level_1_value is not None and current_level_1_value != previous_level_1_value:
                separator_indices.append(i)  # Ajouter l'indice où phylum change

            previous_level_1_value = current_level_1_value

        return separator_indices


def compute_accuracy_from_conf_matrix_df(conf_mat_level):
    correct_predictions = np.diag(conf_mat_level).sum()
    total_predictions = conf_mat_level.sum().sum()
    accuracy = correct_predictions / total_predictions if total_predictions > 0 else 0.0
    return accuracy

