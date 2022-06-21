import matplotlib.pyplot as plt
from json import load
from os import listdir
from statistics import mean
import pandas as pd
import numpy
import seaborn as sns
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd


def load_json(json_file: str) -> dict:
    """
    Charge un fichier json en un dictionnaire
    * json_file (str) : le chemin d'accès au fichier
    """
    return load(open(f"{json_file}", "r"))


def extractor(SPARKLE: str) -> dict[str, str]:
    print("oh.")
    dict_files: dict[str, str] = {}
    for folder in listdir(f"{SPARKLE}/"):
        for subf in listdir(f"{SPARKLE}/{folder}"):
            if 'read' in subf:
                for content in listdir(f"{SPARKLE}/{folder}/{subf}"):
                    if '.json' in content:
                        dict_files = {**dict_files,
                                      folder: f"{SPARKLE}/{folder}/{subf}/{content}"}
    return dict_files


def raw_counts_preds(dict_files, style) -> dict[str, int]:
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    if style == 'one':
        dict_counts: dict[str, int] = {
            **{f"true_{level}": 0 for level in LEVELS}, **{f"false_{level}": 0 for level in LEVELS}, **{f"impossible_{level}": 0 for level in LEVELS}}
    else:
        dict_counts: dict[str, int] = {
            **{f"true_{level}": 0 for level in LEVELS}, **{f"false_{level}": 0 for level in LEVELS}}
    taxa_list: list[str] = list(dict_files.keys())
    for ref, file in dict_files.items():
        temp_storage = load_json(file)
        for i, level in enumerate(LEVELS):
            if ref.split('_')[i] == temp_storage[level]:
                # bien joué ma boi
                dict_counts[f"true_{level}"] = dict_counts[f"true_{level}"] + 1
            elif [x.split('_')[i] for x in taxa_list].count(ref.split('_')[i]) <= 1:
                # has relatives ?
                # dict_counts[f"impossible_{level}"] = dict_counts[f"impossible_{level}"] + 1
                pass
            else:
                # is bad classif :(
                dict_counts[f"false_{level}"] = dict_counts[f"false_{level}"] + 1
    return dict_counts


def good_read_percentage(dict_files) -> dict[str, dict]:
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    taxa_list: list[str] = list(dict_files.keys())
    storage = {k: {} for k in taxa_list}
    for ref, file in dict_files.items():
        temp_storage = load_json(file)
        for i, level in enumerate(LEVELS):
            if f"{ref.split('_')[i]} ({level[0]})" in temp_storage:
                storage[ref][f"{ref.split('_')[i]} ({level[0]})"] = temp_storage[f"{ref.split('_')[i]} ({level[0]})"]
    return storage


def good_read_value(dict_files) -> dict[str, dict]:
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    taxa_list: list[str] = list(dict_files.keys())
    storage = {k: {} for k in taxa_list}
    for ref, file in dict_files.items():
        temp_storage = load_json(file)
        for i, level in enumerate(LEVELS):
            if f"{level} {ref.split('_')[i]}" in temp_storage:
                storage[ref][f"{level} {ref.split('_')[i]}"] = temp_storage[f"{level} {ref.split('_')[i]}"]
    return storage


def bad_read_value(dict_files) -> dict[str, dict]:
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    taxa_list: list[str] = list(dict_files.keys())
    storage = {k: {} for k in taxa_list}
    for ref, file in dict_files.items():
        temp_storage = load_json(file)
        for i, level in enumerate(LEVELS):
            set_all_possibilities = list(set(
                [ref.split('_')[i] for ref in dict_files.keys()]))
            for eltx in set_all_possibilities:
                if f"{level} {eltx}" in temp_storage and eltx != ref.split('_')[i]:
                    storage[ref][f"{level} {ref.split('_')[i]}"] = temp_storage[f"{level} {ref.split('_')[i]}"]
    return storage


def dmatrix_pandas(dict_files):
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    set_all_possibilities = {}
    crosstab_as_list = {}
    for i, level in enumerate(LEVELS):
        set_all_possibilities[level] = list(set(
            [ref.split('_')[i] for ref in dict_files.keys()]))
        crosstab_as_list[level] = [
            [list()]*len(set_all_possibilities[level])]*len(set_all_possibilities[level])
    return set_all_possibilities, crosstab_as_list


def conf_matrix(dict_files):
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    references, tableaux = dmatrix_pandas(dict_files)
    taxa_list: list[str] = list(dict_files.keys())
    for ref, file in dict_files.items():
        temp_storage = load_json(file)
        # gérer ici Tree None (-)
        for i, level in enumerate(LEVELS):
            if f"Tree {ref.split('_')[i]} ({level[0]})" in temp_storage:
                # on a la ref, on la stocke
                for taxa in temp_storage[f"Tree {ref.split('_')[i]} ({level[0]})"]:
                    value = float(temp_storage[taxa][:-1])
                    ref_taxid = references[LEVELS[i+1]
                                           ].index(ref.split('_')[i+1])
                    if taxa[:-4] in references[LEVELS[i+1]]:
                        esm_taxid = references[LEVELS[i+1]].index(taxa[:-4])
                        tableaux[LEVELS[i+1]
                                 ][ref_taxid][esm_taxid].append(value)
    for clade in tableaux.values():
        for col in clade:
            for ligne in col:
                if ligne == []:
                    ligne = [0]
                else:
                    ligne = sum(ligne)/len(ligne)
    print(tableaux['group'])


def good_percentage_aggregator(storage: dict[str, dict]):
    "Rend les pourcentages de bonnes attributions de reads"
    aggregator = {}
    for _, dct in storage.items():
        for key, elt in dct.items():
            if key not in aggregator:
                aggregator[key] = []
            aggregator[key].append(float(elt[:-1]))
    return aggregator


def good_value_aggregator(storage: dict[str, dict]):
    "Rend les pourcentages de bonnes attributions de reads"
    aggregator = {}
    for _, dct in storage.items():
        for key, elt in dct.items():
            if key not in aggregator:
                aggregator[key] = []
            aggregator[key].append(
                float(float(elt.split(' ')[0])/float(elt.split(' ')[4])))
    return aggregator


def read_number_aggregator(storage: dict[str, dict]):
    "Rend les pourcentages de bonnes attributions de reads"
    aggregator_good, aggregator_bad = {}, {}
    for _, dct in storage.items():
        for key, elt in dct.items():
            if key not in aggregator_good:
                aggregator_good[key] = []
            aggregator_good[key].append(int(elt.split(' ')[0]))
            if key not in aggregator_bad:
                aggregator_bad[key] = []
            aggregator_bad[key].append(
                int(elt.split(' ')[4])-int(elt.split(' ')[0]))
    return aggregator_good, aggregator_bad


def int_level(level: str) -> int:
    return ['domain', 'phylum', 'group', 'order', 'family'].index(level)


def ancestor(dim, level, keyset):
    upper_level = int_level(level) - 1
    current_level = int_level(level)
    return [k.split('_')[upper_level] for k in keyset if k.split('_')[current_level] == dim][0]


def crosstab_aggregator(storage: dict[str, str], level: str):
    levels = ['family', 'order', 'group', 'phylum', 'domain']
    # increment de -1 sur le retrieve d'origine
    "Rend les pourcentages de bonnes attributions de reads"
    rast = int_level(level)
    all_sp = list({k.split('_')[rast] for k in storage.keys()})
    df = pd.DataFrame(columns=['id']+all_sp)
    df['id'] = all_sp
    df = df.fillna(0)
    for name, dct in storage.items():
        tpdct = load_json(dct)
        classif_bestiole = name.split('_')[int_level(level)]
        upper_bestiole = name.split('_')[int_level(level)-1]
        if retrieve_parent_level(level) is None:
            list_to_add = [element[:-4]
                           for element in tpdct["Tree None (-)"]]
            for elt in list_to_add:
                if elt in all_sp:
                    df.at[all_sp.index(classif_bestiole), elt] = df.at[all_sp.index(
                        classif_bestiole), elt] + int(tpdct[f"{level} {elt}"].split(' ')[0])
        elif f"Tree {upper_bestiole} ({retrieve_parent_level(level)[0]})" in tpdct:
            list_to_add = [element[:-4]
                           for element in tpdct[f"Tree {upper_bestiole} ({retrieve_parent_level(level)[0]})"]]
            for elt in list_to_add:
                if elt in all_sp:
                    df.at[all_sp.index(classif_bestiole), elt] = df.at[all_sp.index(
                        classif_bestiole), elt] + int(tpdct[f"{level} {elt}"].split(' ')[0])
    df = df.set_index('id')
    ############### P2 : construction de la matrice de clustering réelle ###################
    df_distance_matrix = pd.DataFrame(columns=['id']+all_sp)
    df_distance_matrix['id'] = all_sp
    df_distance_matrix = df_distance_matrix.fillna(1)
    df_distance_matrix = df_distance_matrix.set_index('id')
    for dim_one in all_sp:
        for dim_two in all_sp:
            if dim_one == dim_two:
                df_distance_matrix.at[dim_one,
                                      dim_two] = 0
            ancestor_one = dim_one
            ancestor_two = dim_two
            for enu, lvl in enumerate(levels[levels.index(level):]):
                ancestor_one = ancestor(ancestor_one, lvl, storage.keys())
                ancestor_two = ancestor(ancestor_two, lvl, storage.keys())
                if ancestor_one != ancestor_two:
                    df_distance_matrix.at[dim_one,
                                          dim_two] = df_distance_matrix.at[dim_one, dim_two] + (enu+1)
    return df, df_distance_matrix


def retrieve_parent_level(current: str) -> str | None:
    if int_level(current)-1 >= 0:
        return ['domain', 'phylum', 'group', 'order', 'family'][int_level(current)-1]
    return None


if __name__ == '__main__':
    p = {}
    part = 'all'
    db = 'sampled'
    version = 1
    cols = ['Database', 'true_domain', 'true_phylum', 'true_group', 'true_order', 'true_family', 'false_domain', 'false_phylum', 'false_group', 'false_order', 'false_family', 'impossible_domain', 'impossible_phylum', 'impossible_group',
            'impossible_order', 'impossible_family'] if part == 'one' else ['Database', 'true_domain', 'true_phylum', 'true_group', 'true_order', 'true_family', 'false_domain', 'false_phylum', 'false_group', 'false_order', 'false_family']
    # , f'{part}vsall_sampled_v1', f'{part}vsall_sampled_v2'
    max_number = 1
    for extr in [f'allvsall_{db}_v{version}', f'onevsall_{db}_v{version}']:
        dctx = extractor(extr)
        max_number = len(dctx)
        vals = (good_value_aggregator(good_read_value(dctx)))
        # reads classification accuracy
        averages = {key: mean(val) for key, val in vals.items()}
        # read numbers
        good, bad = (read_number_aggregator(good_read_value(dctx)))
        # print(good)

        tp = []
        for filterd in ['domain', 'phylum', 'group', 'order', 'family']:
            dfd = pd.DataFrame([[k.split(' ')[1]]+[sum(val)]+[sum(bad[k])]
                                for k, val in good.items() if filterd in k], columns=['Taxa', 'Good assignations', 'Bad assignations'])
            dfd = dfd.sort_values('Taxa')
            tp.append(dfd['Good assignations'].sum(
            )/(dfd['Good assignations'].sum()+dfd['Bad assignations'].sum()))
            #print(f"{extr} > {filterd}: {tp[-1]}")
            fig = plt.figure()
            cm = plt.get_cmap('rainbow')
            ax = dfd.plot(x='Taxa',
                          figsize=(15, 6),
                          kind='bar',
                          stacked=True,
                          cmap=cm)
            plt.ylabel('Reads counts')
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2)
            plt.savefig(f"{extr}_{filterd}.png", bbox_inches='tight')
            ##################### seaborn clustermap #######################
            fig, ax = plt.subplots()
            cm = [plt.get_cmap('rainbow'), plt.get_cmap('binary', 4)]
            dendrodf, distMatrix = crosstab_aggregator(dctx, filterd)
            distArray = ssd.squareform(distMatrix)
            distLinkage = hierarchy.linkage(distArray)
            dendrodf = dendrodf.div(dendrodf.sum(axis=1), axis=0) * 100
            for i, y in enumerate([dendrodf]):
                cg = sns.clustermap(y, cmap=cm[i], annot=True, cbar_pos=(.03, 0.1, 0.01, .4),
                                    fmt=".2f", row_linkage=distLinkage, col_linkage=distLinkage)
            # cg.ax_row_dendrogram.set_visible(False)
            #cg2 = sns.clustermap(distMatrix, cmap=cm_b, cbar_pos=None,row_linkage=distLinkage, col_linkage=distLinkage, alpha=0.4)
                ax = cg.ax_heatmap
                ax.set_xlabel('Predicted')
                ax.set_ylabel('Actual')
            plt.ylabel('Percentage of reads')
            plt.savefig(
                f"clustermap_{extr}_{filterd}.png", bbox_inches='tight')

        #print(f"Global : {numpy.prod(tp)}")
        p[''.join(extr.split('_')[0])] = raw_counts_preds(dctx, part)
    df = pd.DataFrame([[kp]+[v for _, v in dt.items()]
                      for kp, dt in p.items()], columns=[f'{db.capitalize()} database (v{version})']+[col.replace('_', ' ') for col in cols[1:]])

    dfp = pd.DataFrame(columns=[f'{db.capitalize()} database (v{version})']+[
                       col.replace('_', ' ') for col in cols[1:]])
    dfp[f'{db.capitalize()} database (v{version})'] = df[f'{db.capitalize()} database (v{version})']
    dfp['true domain'] = (
        df['true domain']/(df['true domain']+df['false domain']))*100
    dfp['false domain'] = (
        df['false domain']/(df['true domain']+df['false domain']))*100
    dfp['true phylum'] = (
        df['true phylum']/(df['true phylum']+df['false phylum']))*100
    dfp['false phylum'] = (
        df['false phylum']/(df['true phylum']+df['false phylum']))*100
    dfp['true group'] = (
        df['true group']/(df['true group']+df['false group']))*100
    dfp['true order'] = (
        df['true order']/(df['true order']+df['false order']))*100
    dfp['false group'] = (
        df['false group']/(df['true group']+df['false group']))*100
    dfp['false order'] = (
        df['false order']/(df['true order']+df['false order']))*100
    dfp['true family'] = (
        df['true family']/(df['true family']+df['false family']))*100
    dfp['false family'] = (
        df['false family']/(df['true family']+df['false family']))*100

    dfp = dfp.round(2)

    fig = plt.figure()
    cm = plt.get_cmap('rainbow', 15)
    ax = dfp.plot(x=f'{db.capitalize()} database (v{version})',
                  figsize=(15, 6),
                  kind='bar',
                  stacked=False,
                  cmap=cm)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim(0, 112)
    plt.ylabel('Percentages of assignations')
    plt.xticks(rotation=0)
    for container in ax.containers:
        labels = [v if v > 0 else "" for v in container.datavalues]
        ax.bar_label(container, labels=labels, rotation=90, padding=3)
    # print(good_read_percentage(dctx))
    """
    print(f"{extr} : {mean(averages.values())}")
    # conf_matrix(dctx)
    plt.bar(range(len(p)), list(p.values()), align='center')
    plt.xticks(range(len(p)), list(p.keys()))

    """
    plt.savefig(f"resume_classif_{part}vsall.png", bbox_inches='tight')
