import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# all_unique_labels = pd.read_csv('all_genome_family.csv',index_col=0)


def plot_conf_mat(ax, conf_mat,  level, filename=None):

    # Déterminer si on doit afficher les labels
    afficher_labels = len(conf_mat) <= 15
    # Création de la figure
    cmap = plt.cm.RdBu_r  # Colormap avec bon contraste
    norm = mcolors.PowerNorm(gamma=0.5)  # Accentue les différences des faibles valeurs
    common_names = [taxon for taxon in all_unique_labels['family'] if taxon in conf_mat_level.index]
    sub_conf_mat = conf_mat_level.loc[common_names, common_names]

    title = f'Conf_mat for level {level}'

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
    # else:
    #     plt.show()

if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(8, 6))
    level = 'group'
    csv_file = f"/home/hcourtei/Projects/MicroTaxo/codes/genouest_archive/model_base_02_05_11_26/eval/metrics/ConfMat_{level}.csv"

    conf_mat_level = pd.read_csv(csv_file, index_col=0, sep= ";")

    print(conf_mat_level)
    print(all_unique_labels['group'].unique())
    # plot_conf_mat(ax, conf_mat_level,level=level, filename= None) #csv_file.replace(".csv", ""))
    # plt.show()
    common_names = [taxon for taxon in all_unique_labels[level] if taxon in conf_mat_level.index]
    sub_conf_mat = conf_mat_level.loc[common_names, common_names]
    # print(sub_conf_mat)
    # plot_conf_mat(ax, sub_conf_mat,level=level, filename= None) #csv_file.replace(".csv", ""))

    # #
    # # # Par exemple, si votre mapping provient de votre CSV :
    # # mapping_phylum = all_unique_labels.set_index('family')['phylum'].to_dict()
    # #
    # # # Extraire la liste des phyla dans l'ordre des labels de la matrice
    # # phylum_list = [mapping_phylum[label] for label in sub_conf_mat.index if label in mapping_phylum]
    # #
    # # # Identifier les indices où le phylum change
    # # sep_indices = [i for i in range(1, len(phylum_list)) if phylum_list[i] != phylum_list[i - 1]]
    # #
    # # for pos in sep_indices:
    # #     ax.axhline(pos, color='black', linewidth=2)
    # #     ax.axvline(pos, color='black', linewidth=2)
    # # plt.show()
    # # Créer le mapping family -> phylum et normaliser les chaînes si nécessaire
    # mapping_phylum = all_unique_labels.set_index('family')['phylum'].to_dict()
    # mapping_phylum = {k.strip(): v.strip() for k, v in mapping_phylum.items()}
    # #
    # # Filtrer les familles présentes dans la matrice
    # common_names = [taxon for taxon in all_unique_labels[level] if taxon in conf_mat_level.index]
    # sub_conf_mat = conf_mat_level.loc[common_names, common_names]
    #
    # # Trier la sous-matrice par phylum pour que les familles du même phylum soient groupées
    # sub_conf_mat = sub_conf_mat.sort_index(key=lambda x: [mapping_phylum.get(i.strip(), "") for i in x])
    #
    # # Extraire la liste des phyla dans l'ordre des labels
    # phylum_list = [mapping_phylum.get(label.strip(), "") for label in sub_conf_mat.index]
    # print("Phyla dans l'ordre :", phylum_list)
    # print("Phyla uniques :", list(dict.fromkeys(phylum_list)))
    #
    # # Identifier les indices où le phylum change (ce qui correspond aux frontières entre groupes)
    # sep_indices = []
    # prev = phylum_list[0]
    # for i, ph in enumerate(phylum_list[1:], start=1):
    #     if ph != prev:
    #         sep_indices.append(i)
    #         prev = ph
    # print("Indices de séparation :", sep_indices)
    # for pos in sep_indices:
    #     ax.axhline(pos, color='black', linewidth=2)
    #     ax.axvline(pos, color='black', linewidth=2)
    # plt.show()

    # for level in ["phylum", "group", "order", "family" ]:
    #     print("-"*30)
    #     print("duplicates for level : ", level)
    #     l= list(all_unique_labels[level][all_unique_labels[level].duplicated()].unique())
    #     print(l)

    orders_paths = all_unique_labels[['domain', 'phylum', 'group', 'order']].drop_duplicates()
    duplicated_orders = orders_paths.groupby('order').filter(lambda x: len(x) > 1)
    duplicated_orders.sort_values('order')
    group_paths = all_unique_labels[['domain', 'phylum', 'group']].drop_duplicates()
    duplicated_groups = group_paths.groupby('group').filter(lambda x: len(x) > 1)
    duplicated_groups.sort_values('group')
    phylum_paths = all_unique_labels[['domain', 'phylum']].drop_duplicates()
    phylum_groups = phylum_paths.groupby('phylum').filter(lambda x: len(x) > 1)
    phylum_groups.sort_values('phylum')

