import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_conf_mat(conf_mat,  title=f"Conf Mat for level xxx", filename=None):

    # Déterminer si on doit afficher les labels
    afficher_labels = len(conf_mat) <= 15
    # Création de la figure
    cmap = plt.cm.RdBu_r  # Colormap avec bon contraste
    norm = mcolors.PowerNorm(gamma=0.5)  # Accentue les différences des faibles valeurs

    fig, ax = plt.subplots(figsize=(8, 6))
    if afficher_labels:
        sns.heatmap(conf_mat, annot=True, fmt="d", cmap="Blues",
                    xticklabels=afficher_labels, yticklabels=afficher_labels, ax=ax)
        plt.xticks(rotation=45, ha="right")  # Rotation des prédictions
        plt.yticks(rotation=45, va="top")  # Rotation des réels
    else:
        cax = ax.imshow(conf_mat, cmap=cmap, norm=norm, interpolation="nearest")    # Ajout des titres
        plt.colorbar(cax, ax=ax)

    plt.title(title)
    plt.xlabel("Prédit")
    plt.ylabel("Réel")
    plt.tight_layout()
    # Affichage de la matrice
    if filename:
        plt.savefig(filename, dpi=300)
        # Affichage de la matrice
    else:
        plt.show()

if __name__ == '__main__':

    csv_file = "/home/hcourtei/Projects/MicroTaxo/codes/genouest_archive/model_base_02_04_15_44/metrics/ConfMat_phylum.csv"
    level = 'phylum'
    conf_mat_level = pd.read_csv(csv_file, index_col=0, sep= ";")
    print(conf_mat_level)
    plot_conf_mat(conf_mat_level, title=f"Conf Mat for level {level}", filename=csv_file.replace(".csv", ""))
