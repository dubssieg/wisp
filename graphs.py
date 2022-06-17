import matplotlib.pyplot as plt
from json import load
from os import listdir
from statistics import mean
import pandas as pd


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
                dict_counts[f"impossible_{level}"] = dict_counts[f"impossible_{level}"] + 1
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


if __name__ == '__main__':
    p = {}
    part = 'all'
    cols = ['Database', 'true_domain', 'true_phylum', 'true_group', 'true_order', 'true_family', 'false_domain', 'false_phylum', 'false_group', 'false_order', 'false_family', 'impossible_domain', 'impossible_phylum', 'impossible_group',
            'impossible_order', 'impossible_family'] if part == 'one' else ['Database', 'true_domain', 'true_phylum', 'true_group', 'true_order', 'true_family', 'false_domain', 'false_phylum', 'false_group', 'false_order', 'false_family']
    # , f'{part}vsall_subsampled', f'{part}vsall_subsampled_v2'
    for extr in [f'{part}vsall_small_v1', f'{part}vsall_small_v2', f'{part}vsall_sampled_v1', f'{part}vsall_sampled_v2']:
        dctx = extractor(extr)
        vals = (good_value_aggregator(good_read_value(dctx)))
        # reads classification accuracy
        averages = {key: mean(val) for key, val in vals.items()}
        # read numbers
        good, bad = (read_number_aggregator(good_read_value(dctx)))
        for filterd in ['domain', 'phylum', 'group', 'order', 'family']:
            dfd = pd.DataFrame([[k.split(' ')[1]]+[sum(val)]+[sum(bad[k])]
                                for k, val in good.items() if filterd in k], columns=['Taxa', 'Good assignations', 'Bad assignations'])
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

        p[' '.join(extr.split('_')[1:])] = raw_counts_preds(dctx, part)
    df = pd.DataFrame([[kp]+[v for _, v in dt.items()]
                      for kp, dt in p.items()], columns=['Database']+[col.replace('_', ' ') for col in cols[1:]])
    fig = plt.figure()
    cm = plt.get_cmap('rainbow', 15)
    ax = df.plot(x='Database',
                 figsize=(15, 6),
                 kind='bar',
                 stacked=False,
                 cmap=cm)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim(0, 125)
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
