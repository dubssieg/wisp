import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_conf_mat(conf_mat,  title=f"Conf Mat for level xxx", filename=None):

    # Déterminer si on doit afficher les labels
    afficher_labels = len(conf_mat) <= 20
    # Création de la figure
    plt.figure(figsize=(8, 6))
    sns.heatmap(conf_mat, annot=True, fmt="d", cmap="Blues",
                xticklabels=afficher_labels,
                yticklabels=afficher_labels)

    # Ajout des titres
    plt.title(title)
    plt.xlabel("Prédit")
    plt.ylabel("Réel")
    if afficher_labels:
        plt.xticks(rotation=45, ha="right")  # Rotation des prédictions
        plt.yticks(rotation=45, va="top")  # Rotation des réels

    plt.tight_layout()
    # Affichage de la matrice
    if filename:
        plt.savefig(filename, dpi=300)
        # Affichage de la matrice
    else:
        plt.show()
