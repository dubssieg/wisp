from collections import defaultdict
from sklearn.metrics import confusion_matrix
import numpy as np
import pandas as pd


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
        for level in ['domain', 'phylum', 'group', 'order', 'family']:
            true_value = true_labels.get(level)
            pred_value = pred_labels.get(level)

            if true_value is not None and pred_value is not None:  # Exclure les None
                self.true_labels[level].append(true_value)
                self.pred_labels[level].append(pred_value)

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

        y_true = self.true_labels[level]
        y_pred = self.pred_labels[level]

        # Filtrer les None
        filtered_true = [label for label in y_true if label is not None]
        filtered_pred = [label for label in y_pred if label is not None]

        # Liste des classes uniques sans None
        classes = sorted(set(filtered_true + filtered_pred))

        if not classes:
            return pd.DataFrame()  # Matrice vide si aucun taxon valide

        # matrix = confusion_matrix(filtered_true, filtered_pred, labels=classes) # replace with
        all_classes = sorted(set(self.true_labels[level]) | set(self.pred_labels[level]))

        # Calculer la matrice de confusion avec toutes les classes connues
        matrix = confusion_matrix(filtered_true, filtered_pred, labels=all_classes)

        return pd.DataFrame(matrix, index=classes, columns=all_classes) # columns=all_classes)

    def get_all_confusion_matrices(self):
        """
        Retourne un dictionnaire contenant les matrices de confusion pour tous les niveaux taxonomiques.

        Returns:
            dict: {niveau: matrice de confusion sous forme de DataFrame}
        """
        return {level: self.get_confusion_matrix(level) for level in self.true_labels if self.true_labels[level]}


def compute_accuracy_from_conf_matrix_df(conf_mat_level):
    # Somme des éléments de la diagonale (prédictions correctes)
    correct_predictions = np.diag(conf_mat_level).sum()

    # Somme de tous les éléments de la matrice (total des prédictions)
    total_predictions = conf_mat_level.sum().sum()

    # Calcul de l'accuracy
    accuracy = correct_predictions / total_predictions if total_predictions > 0 else 0.0

    return accuracy