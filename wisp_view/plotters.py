"Analysis upon kmers repartition and stuff"

from argparse import ArgumentParser
from string import ascii_lowercase
import matplotlib.pyplot as plt
import pandas as pd
from wisp_lib import kmer_indexing_brut, recode_kmer_4, kmer_indexing_10000
from wisp_tools import my_parser
from collections import Counter
import xgboost as xgb
from os import listdir
from statistics import mean, stdev
import numpy as np
import matplotlib.pyplot as plt
from json import load
from os import listdir
from statistics import mean
import pandas as pd
import numpy
import seaborn as sns
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd
import pylab as pl
import numpy as np
import matplotlib as mpl
from pathlib import Path
from mpl_toolkits import mplot3d
from itertools import chain


def load_json(json_file: str) -> dict:
    """
    Charge un fichier json en un dictionnaire
    * json_file (str) : le chemin d'accès au fichier
    """
    return load(open(f"{json_file}", "r"))


def extractor(SPARKLE: str) -> dict[str, str]:
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


def cbar_plotting():
    a = np.array([[0, 100]])
    pl.figure(figsize=(9, 1.5))
    img = pl.imshow(a, cmap="rainbow")
    pl.gca().set_visible(False)
    cax = pl.axes([0.1, 0.2, 0.8, 0.6])
    pl.colorbar(orientation="horizontal", cax=cax)
    pl.savefig("output/figures/clustering/colorbar.png")


def cbar_plotting_2():
    fig, ax = plt.subplots(1, 1)
    fraction = 1  # .05
    norm = mpl.colors.Normalize(vmin=0, vmax=100)
    cbar = ax.figure.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap='rainbow'),
        ax=ax, pad=.05, extend='both', fraction=fraction)
    ax.axis('off')
    plt.savefig("output/figures/clustering/cbar.png", bbox_inches='tight')


def clustering_plotting(folders_to_analyze: list[str]):
    plt.rcParams.update({'figure.max_open_warning': 0})
    cbar_plotting()
    cbar_plotting_2()
    p = {}
    dbs = ['supported']  # , 'supported'
    for db in dbs:
        version = 1
        cols = ['Database', 'true_domain', 'true_phylum', 'true_group', 'true_order', 'true_family',
                'false_domain', 'false_phylum', 'false_group', 'false_order', 'false_family']
        #
        max_number = 1
        # [f'baseline_{db}_v{version}', f'leaveoneout_{db}_v{version}', f'baseline6percent_{db}_v{version}', f'leaveoneout6percent_{db}_v{version}']:
        for extr in [f'leaveoneout_bernard_v{version}', f'leaveoneout6percent_bernard_v{version}', f'leaveoneout_supported_v{version}', f'leaveoneout6percent_supported_v{version}']:
            dctx = extractor(extr)
            max_number = len(dctx)
            vals = (good_value_aggregator(good_read_value(dctx)))
            # reads classification accuracy
            averages = {key: mean(val) for key, val in vals.items()}
            # read numbers
            good, bad = (read_number_aggregator(good_read_value(dctx)))

            tp = []
            for filterd in ['domain', 'phylum', 'group', 'order', 'family']:
                Sn_list, Sp_list = [], []
                dfd = pd.DataFrame([[k.split(' ')[1]]+[sum(val)]+[sum(bad[k])]
                                    for k, val in good.items() if filterd in k], columns=['Taxa', 'Good assignations', 'Bad assignations'])
                dfd = dfd.sort_values('Taxa')
                tp.append(dfd['Good assignations'].sum(
                )/(dfd['Good assignations'].sum()+dfd['Bad assignations'].sum()))
                print("")
                print(f"{extr} ({filterd})")
                print(f"Accuracy : {tp[-1]}")
                fig = plt.figure()
                cm = plt.get_cmap('rainbow')
                ax = dfd.plot(x='Taxa',
                              figsize=(15, 6),
                              kind='bar',
                              stacked=True,
                              cmap=cm)
                plt.ylabel('Reads counts')
                plt.legend(loc='upper center',
                           bbox_to_anchor=(0.5, 1.1), ncol=2)
                plt.savefig(
                    f"output/figures/clustering/{extr}_{filterd}.png", bbox_inches='tight')
                ##################### seaborn clustermap #######################
                fig, ax = plt.subplots()
                cm = [plt.get_cmap('rainbow'), plt.get_cmap('binary', 4)]
                dendrodf, distMatrix = crosstab_aggregator(dctx, filterd)
                dendrodf = dendrodf.transpose()
                classes = list(dendrodf.index)
                for i, classe in enumerate(classes):
                    print(
                        f"{classe} : Sn = {dendrodf.iloc[i][classe]/sum(dendrodf.iloc[i])} ; Sp = {(sum(np.diag(dendrodf))-dendrodf.iloc[i][classe])/((sum(np.diag(dendrodf))-dendrodf.iloc[i][classe])+(sum(dendrodf[classe])-dendrodf.iloc[i][classe]))}")
                    Sn_list.append(
                        dendrodf.iloc[i][classe]/sum(dendrodf.iloc[i]))
                    Sp_list.append((sum(np.diag(dendrodf))-dendrodf.iloc[i][classe])/((sum(np.diag(
                        dendrodf))-dendrodf.iloc[i][classe])+(sum(dendrodf[classe])-dendrodf.iloc[i][classe])))
                print(f"Average Sn ({filterd}) : {mean(Sn_list)}")
                print(f"Average Sp ({filterd}) : {mean(Sp_list)}")
                distArray = ssd.squareform(distMatrix)
                distLinkage = hierarchy.linkage(distArray)
                dendrodf = dendrodf.div(dendrodf.sum(axis=1), axis=0) * 100
                for i, y in enumerate([dendrodf]):
                    cg = sns.clustermap(y, cmap=cm[i], cbar_pos=None, mask=distMatrix > 1,
                                        row_linkage=distLinkage, col_linkage=distLinkage, xticklabels=True, yticklabels=True, annot=True, fmt=".2f")  #
                # (1, 0.2, 0.01, .6)
                # cg.ax_row_dendrogram.set_visible(False)
                # cg2 = sns.clustermap(distMatrix, cmap=cm_b, cbar_pos=None,row_linkage=distLinkage, col_linkage=distLinkage, alpha=0.4)
                    ax = cg.ax_heatmap
                    ax.set_facecolor("#d6dbdf")
                    ax.tick_params(labelsize=7)
                    ax.set_xlabel('Predicted')
                    ax.set_ylabel('Actual')

                    for t in ax.texts:
                        if float(t.get_text()) > 0.0:
                            t.set_text(t.get_text())
                        else:
                            t.set_text("")
                plt.ylabel('Actual')  # Percentage of reads
                plt.savefig(
                    f"output/figures/clustering/clustermap_{extr}_{filterd}_percentages.png", bbox_inches='tight')
                ##################### seaborn clustermap #######################
                fig, ax = plt.subplots()
                cm = [plt.get_cmap('rainbow')]
                for i, y in enumerate([dendrodf]):
                    cg = sns.clustermap(y, cmap=cm[0], cbar_pos=None, mask=distMatrix > 1,
                                        row_linkage=distLinkage, col_linkage=distLinkage, xticklabels=True, yticklabels=True)  # , annot=True, fmt=".2f"
                # (1, 0.2, 0.01, .6)
                # cg.ax_row_dendrogram.set_visible(False)
                # cg2 = sns.clustermap(distMatrix, cmap=cm_b, cbar_pos=None,row_linkage=distLinkage, col_linkage=distLinkage, alpha=0.4)
                    ax = cg.ax_heatmap
                    ax.set_facecolor("#d6dbdf")
                    ax.tick_params(labelsize=7)
                    ax.set_xlabel('Predicted')
                    ax.set_ylabel('Actual')

                    for t in ax.texts:
                        if float(t.get_text()) > 0.0:
                            t.set_text(t.get_text())
                        else:
                            t.set_text("")
                plt.ylabel('Actual')  # Percentage of reads
                plt.savefig(
                    f"clustermap_{extr}_{filterd}_nopercentages_transparent.png", bbox_inches='tight', transparent=True)

            print(f"Cumulative accuracy : {numpy.prod(tp)}")
            p[''.join(extr.split('_')[0])+" "+''.join(extr.split('_')
                                                      [1])] = raw_counts_preds(dctx, 'all')
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

        print(dfp)

        fig = plt.figure()
        cm = plt.get_cmap('rainbow', 15)
        ax = dfp.plot(figsize=(15, 6),
                      kind='bar',
                      stacked=False,
                      cmap=cm)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim(0, 112)
        plt.ylabel('Percentages of assignations')
        plt.xticks([0, 1, 2, 3], ['Bernard $v1$', 'Bernard $v1E$',
                   'Supported $v1$', 'Supported $v1E$'], rotation=0)

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
        plt.savefig("leaveoneout_global_results.png", bbox_inches='tight')


def my_encoder_4():
    # maybe this can help gain speed ? specific to k=4
    return {f"{a}{b}{c}{d}": 0 for a in ['A', 'T', 'G', 'C'] for b in ['A', 'T', 'G', 'C'] for c in ['A', 'T', 'G', 'C']for d in ['A', 'T', 'G', 'C']}


def signatures(list_sequences):
    alphabet = Counter(my_encoder_4())
    total_len = 0
    for sequence in list_sequences:
        total_len += len(sequence)
        alphabet += kmer_indexing_brut(sequence, 4)
    return Counter({key: float(value/total_len) for key, value in alphabet.items()})


def signatures_sd(list_sequences):
    alphabet = {k: [] for k in my_encoder_4().keys()}
    for sequence in list_sequences:
        for kmer, kval in kmer_indexing_10000(sequence, 4).items():
            alphabet[kmer].append(kval)
    return Counter({key: float(mean(value)) if len(value) > 1 else value for key, value in alphabet.items()})


def code(value, mn, sd):
    if value < mn - sd:
        return 0
    elif value > mn + sd:
        return 2
    else:
        return 1


def compute_signatures(level, pwd, listing):
    rets, raw_rets, dev_rets = {}, {}, {}
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    splitting = LEVELS.index(level)
    targets = list(set(s.split('_')[splitting] for s in listdir(pwd)))
    for t in targets:
        all_genomes_matching = [
            f"{pwd}/{tr}" for tr in listdir(pwd) if tr.split('_')[splitting] == t]
        list_sequences = [my_parser(genome_matching, True, True, 'm')[
            'm'] for genome_matching in all_genomes_matching]
        res = signatures(list_sequences)
        deviation = signatures_sd(list_sequences)
        mn = mean(res.values())
        sd = stdev(res.values())
        tpd = {k: code(v, mn, sd) for k, v in res.items()}
        tabl, tabl2, tabldev = [['' for _ in range(16)] for _ in range(16)], [
            ['' for _ in range(16)] for _ in range(16)], [
            ['' for _ in range(16)] for _ in range(16)]
        for k, v in tpd.items():
            tabl[listing.index(k[:2])][listing.index(k[2:])] = v
        for k, v in res.items():
            tabl2[listing.index(k[:2])][listing.index(k[2:])] = v
        for k, v in deviation.items():
            tabldev[listing.index(k[:2])][listing.index(k[2:])] = v
        rets[t] = np.asarray(tabl)
        raw_rets[t] = np.asarray(tabl2)
        dev_rets[t] = np.asarray(tabldev)

    return rets, raw_rets, dev_rets


def plot_database_features(db_path, output_path):  # ex : data/small/
    listing = [f"{a}{b}" for a in ['A', 'T', 'G', 'C']
               for b in ['A', 'T', 'G', 'C']]
    Path(f"{output_path}/").mkdir(parents=True, exist_ok=True)
    for level in listdir(db_path):
        files = []
        if level in ['family', 'group', 'order']:
            Path(f"{output_path}/{level}").mkdir(parents=True, exist_ok=True)
            Path(
                f"{output_path}/{level}/features_3d").mkdir(parents=True, exist_ok=True)
            files = [f"{db_path}/{level}/{file}" for file in listdir(
                f"{db_path}/{level}") if 'saved_model.json' in file]
            for file in files:
                print(f"Making {file.split('/')[-1].split('_')[0]} graph!")
                plot_some_features(file, listing, file.split(
                    '/')[-1].split('_')[0], f"{output_path}/{level}")


def plot_some_features(my_path, listing, filename, output_path):
    bst = xgb.Booster()
    bst.load_model(my_path)
    mapped = {recode_kmer_4(str(k[1:]), 4): v for k, v in bst.get_score(
        importance_type='gain').items()}

    # then we plot features in 3d
    sds = [[0 for _ in range(16)] for _ in range(16)]
    for k, v in mapped.items():
        sds[listing.index(k[:2])][listing.index(k[2:])] = float(v)
    plt.figure()
    cm = plt.get_cmap('rainbow')
    ax3 = plt.axes(projection='3d')
    x = np.arange(0, 16, 1)
    y = np.arange(0, 16, 1)
    X, Y = np.meshgrid(x, y)
    for maskd in ['AA', 'TT', 'CC', 'GG']:
        idx = listing.index(maskd)
        sds[idx][idx] = np.nan
    sds = np.asarray(sds)
    ax3.set_title(my_path.split('/')[-1].split('_')[0])
    ax3.set_xticks([i for i in range(16)])
    ax3.set_yticks([i for i in range(16)])
    ax3.set_xticklabels([listing[i] for i in range(16)])
    ax3.set_yticklabels([listing[i] for i in range(16)])
    ax3.view_init(60, 35)
    ax3.tick_params(axis=u'both', which=u'both', length=0)
    ax3.plot_surface(X, Y, sds, cmap=cm, edgecolor='none')
    plt.savefig(f"{output_path}/features_3d/{filename}_3d_features.png",
                bbox_inches='tight', transparent=True)


def plot_repartition_top_kmers(number_to_plot: int, sequence: str, ksize: int, output_path: str) -> None:
    # gives most common at global scale
    counter = kmer_indexing_brut(
        sequence, ksize).most_common(number_to_plot)
    elements = [e for (e, _) in counter]
    df = pd.DataFrame(columns=elements)
    # for those, we compute locally their abundance
    for i in range(0, len(sequence), int(len(sequence)/50000)):
        local_kmers = kmer_indexing_brut(sequence[i:i+50000], ksize)
        my_local_kmers = {k: v for k,
                          v in local_kmers.items() if k in elements}
        label = f"[{i}:{i+50000}]"
        my_series = pd.Series(data=my_local_kmers, index=list(
            my_local_kmers.keys()), name=label)
        df = df.append(my_series).fillna(0)
    df = df.transpose()
    ax = df.plot(figsize=(
        20, 6), ylabel=f'kmers/{int(len(sequence)/5000)}bp', rot=90, colormap='cividis')
    plt.savefig(f"{output_path}/kmers_repartition.svg", bbox_inches='tight')


def delta_sequence(seq1: str, seq2: str, ksize: int, output_path: str) -> None:
    counts_1, counts_2 = kmer_indexing_brut(
        seq1, ksize), kmer_indexing_brut(seq2, ksize)
    counts_1, counts_2 = Counter({k: v/len(seq1) for k, v in counts_1.items()}), Counter({
        k: v/len(seq2) for k, v in counts_2.items()})
    counts_1.subtract(counts_2)
    figure = plt.figure(figsize=(7, 20))
    plt.plot(counts_1.values(), counts_1.keys())
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.savefig(f"{output_path}/delta_kmers.png", bbox_inches='tight')


def list_retirever(key, d: list[dict]) -> list:
    my_list = []
    for di in d:
        if key in di:
            my_list += [di[key]]
        else:
            my_list += [0]
    return my_list


def plot_stacked_values(di: dict[str, dict], dref: dict, output_path: str) -> None:
    temp_counter = sum((list(di.values())[0]).values())
    # for k, v in dref.items(): dref[k] = int((v*temp_counter)/100)
    di['Reference'] = dref
    xi = list(di.keys())
    aggregate_keys = []
    for d in di.values():
        aggregate_keys += d.keys()
    aggregate_keys = set(aggregate_keys)
    sample_keys = set(dref.keys())
    cm = sns.color_palette('rainbow', len(sample_keys)+5, as_cmap=False)
    cm_metagenomic = sns.color_palette(
        'flare', len(aggregate_keys), as_cmap=False)

    for key, rawdict in di.items():
        if key != 'Reference':
            di[key] = {k: float(v/sum(rawdict.values())*100)
                       for k, v in rawdict.items()}

    aggregator_falsepositive = {k: np.array(list_retirever(k, list(di.values())))
                                for k in aggregate_keys-sample_keys}
    aggregator_truepositive = {k: np.array(list_retirever(k, list(di.values())))
                               for k in sample_keys}
    all_others_y = np.array([0.0 for _ in range(len(di))])

    figure = plt.figure(figsize=(2*len(di), 25))
    for i, (label, y) in enumerate(aggregator_truepositive.items()):
        plt.bar(xi, y, label=label, bottom=all_others_y, color=cm[i])
        all_others_y += y
    for i, (label, y) in enumerate(aggregator_falsepositive.items()):
        plt.bar(xi, y, label=label, bottom=all_others_y,
                color=cm_metagenomic[i])
        all_others_y += y
    plt.bar(xi, np.array([0 for _ in range(len(di)-1)]+[100-sum(di['Reference'].values())]),
            label="Without barcode", bottom=all_others_y, color='grey')
    plt.legend(loc='center left', bbox_to_anchor=(
        1, 0.5), frameon=False)
    plt.yticks(rotation=270)
    plt.savefig(f"{output_path}/compare_outputs.png", bbox_inches='tight')


def mfunc(x):
    try:
        return stdev(x)
    except TypeError:
        return x


def compdiff_plotting(input_dir, output_path):
    listing = [f"{a}{b}" for a in ['A', 'T', 'G', 'C']
               for b in ['A', 'T', 'G', 'C']]
    Path(f"{output_path}/").mkdir(parents=True, exist_ok=True)
    for level in ['domain', 'phylum', 'group', 'order', 'family']:
        Path(f"{output_path}/{level}").mkdir(parents=True, exist_ok=True)
        Path(f"{output_path}/{level}/sigmadiff").mkdir(parents=True, exist_ok=True)
        Path(f"{output_path}/{level}/compdiff").mkdir(parents=True, exist_ok=True)
        Path(f"{output_path}/{level}/compdiff_3d").mkdir(parents=True, exist_ok=True)
        elts, raw_elts, dev_items = compute_signatures(
            level, input_dir, listing)
        for key, elt in dev_items.items():
            fig = plt.figure()
            cm = plt.get_cmap('rainbow')
            ax = fig.add_subplot(111)
            cax = ax.matshow(elt, cmap=cm)
            # plt.style.use('dark_background') vmin=0, vmax=2
            plt.title(f"{key}")
            plt.xticks(rotation=90)
            plt.yticks(fontsize=9)
            plt.xticks(fontsize=9)
            ax.set_xticks([i for i in range(16)])
            ax.set_yticks([i for i in range(16)])
            ax.set_xticklabels([listing[i] for i in range(16)])
            ax.set_yticklabels([listing[i] for i in range(16)])
            ax.tick_params(axis=u'both', which=u'both', length=0)
            cbar = fig.colorbar(cax)
            cbar.ax.set_title(
                '$\sigma_{freq}$')
            # cbar.ax.set_yticks([0, 1, 2])
            # cbar.ax.set_yticklabels(['$f < \mu - \sigma$', '$f = \mu \pm \sigma$', '$f > \mu + \sigma$'])
            plt.savefig(f"{output_path}/{level}/sigmadiff/{key}_sigmadiff.png",
                        bbox_inches='tight', transparent=True)
        for key, elt in elts.items():
            print("\n"+key+"\n")
            print(elt)
            fig = plt.figure()
            cm = plt.get_cmap('rainbow', 3)
            ax = fig.add_subplot(111)
            cax = ax.matshow(elt, cmap=cm, vmin=0, vmax=2)
            # plt.style.use('dark_background')
            plt.title(f"{key}")
            plt.xticks(rotation=90)
            plt.yticks(fontsize=9)
            plt.xticks(fontsize=9)
            ax.set_xticks([i for i in range(16)])
            ax.set_yticks([i for i in range(16)])
            ax.set_xticklabels([listing[i] for i in range(16)])
            ax.set_yticklabels([listing[i] for i in range(16)])
            ax.tick_params(axis=u'both', which=u'both', length=0)
            cbar = fig.colorbar(cax)
            cbar.ax.set_yticks([0, 1, 2])
            cbar.ax.set_yticklabels(
                ['$f < \mu - \sigma$', '$f = \mu \pm \sigma$', '$f > \mu + \sigma$'])
            plt.savefig(f"{output_path}/{level}/compdiff/{key}_compdiff.png",
                        bbox_inches='tight', transparent=True)
        for key, elt in raw_elts.items():
            fig3 = plt.figure()
            cm = plt.get_cmap('rainbow')
            ax3 = plt.axes(projection='3d')
            x = np.arange(0, 16, 1)
            y = np.arange(0, 16, 1)
            X, Y = np.meshgrid(x, y)
            sds = [['' for _ in range(16)] for _ in range(16)]
            for i, x_axis in enumerate(elt):
                for j, y_axis in enumerate(x_axis):
                    sds[i][j] = float(y_axis)
            for maskd in ['AA', 'TT', 'CC', 'GG']:
                idx = listing.index(maskd)
                sds[idx][idx] = np.nan
            sds = np.asarray(sds)
            ax3.set_title(key)
            ax3.set_xticks([i for i in range(16)])
            ax3.set_yticks([i for i in range(16)])
            ax3.set_xticklabels([listing[i] for i in range(16)])
            ax3.set_yticklabels([listing[i] for i in range(16)])
            ax3.view_init(60, 35)
            ax3.tick_params(axis=u'both', which=u'both', length=0)
            ax3.plot_surface(X, Y, sds, cmap=cm, edgecolor='none')
            plt.savefig(f"{output_path}/{level}/compdiff_3d/{key}_3d_compdiff.png",
                        bbox_inches='tight', transparent=True)
        eltx = np.dstack(list(raw_elts.values()))
        fig2 = plt.figure()
        cm = plt.get_cmap('rainbow')
        ax2 = plt.axes(projection='3d')
        x = np.arange(0, 16, 1)
        y = np.arange(0, 16, 1)
        sds = [['' for _ in range(16)] for _ in range(16)]
        for i, x_axis in enumerate(eltx):
            for j, y_axis in enumerate(x_axis):
                sds[i][j] = stdev(y_axis) if len(y_axis) > 1 else float(y_axis)
        for maskd in ['AA', 'TT', 'CC', 'GG']:
            idx = listing.index(maskd)
            sds[idx][idx] = np.nan
        sds = np.asarray(sds)
        # np.vectorize(mfunc)(eltx)
        X, Y = np.meshgrid(x, y)
        ax2.set_xticks([i for i in range(16)])
        ax2.set_yticks([i for i in range(16)])
        ax2.set_xticklabels([listing[i] for i in range(16)])
        ax2.set_yticklabels([listing[i] for i in range(16)])
        ax2.set_title("All genomes")
        ax2.view_init(60, 35)
        ax2.tick_params(axis=u'both', which=u'both', length=0)
        ax2.plot_surface(X, Y, sds, cmap=cm, edgecolor='none')
        plt.savefig(f"{output_path}/database_3d_compdiff.png",
                    bbox_inches='tight', transparent=True)
