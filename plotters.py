"Analysis upon kmers repartition and stuff"

import matplotlib.pyplot as plt
import pandas as pd
from wisp_lib import kmer_indexing_brut, recode_kmer_4
from wisp_view import plot_features
from python_tools import my_parser
from collections import Counter
import xgboost as xgb
from os import listdir
from statistics import mean, stdev
import numpy as np


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


def code(value, mn, sd):
    if value < mn - sd:
        return 0
    elif value > mn + sd:
        return 2
    else:
        return 1


def compute_signatures(level, pwd, listing):
    rets = {}
    LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    splitting = LEVELS.index(level)
    targets = list(set(s.split('_')[splitting] for s in listdir(pwd)))
    for t in targets:
        all_genomes_matching = [
            f"{pwd}/{tr}" for tr in listdir(pwd) if tr.split('_')[splitting] == t]
        list_sequences = [my_parser(genome_matching, True, True, 'm')[
            'm'] for genome_matching in all_genomes_matching]
        res = signatures(list_sequences)
        mn = mean(res.values())
        sd = stdev(res.values())
        tpd = {k: code(v, mn, sd) for k, v in res.items()}
        tabl = [['' for _ in range(16)] for _ in range(16)]
        for k, v in tpd.items():
            tabl[listing.index(k[:2])][listing.index(k[2:])] = v
        rets[t] = np.asarray(tabl)
    return rets


def plot_database_features(db_path):  # ex : data/small/
    files = []
    db_path = f"data/{db_path}/"
    for level in listdir(db_path):
        files.extend([f"{db_path}/{level}/{file}" for file in listdir(
            f"{db_path}/{level}") if 'saved_model.json' in file])
    for file in files:
        plot_some_features(file)


def plot_some_features(my_path):
    bst = xgb.Booster()
    bst.load_model(my_path)
    mapped = {
        f"{recode_kmer_4(str(k[1:]),max([len(str(ky[1:])) for ky, _ in bst.get_score(importance_type='gain').items()]))}": v for k, v in bst.get_score(importance_type='gain').items()}
    if mapped != {}:
        plot_features(f"{OUTPUT_PATH}{my_path.split('/')[-1].split('_')[0]}_", mapped, "",
                      "", "")


def plot_repartition_top_kmers(number_to_plot: int, sequence: str, ksize: int) -> None:
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
    plt.savefig(f"{OUTPUT_PATH}kmers_repartition.svg", bbox_inches='tight')


def delta_sequence(seq1: str, seq2: str, pattern: str, ksize: int) -> None:
    counts_1, counts_2 = kmer_indexing_brut(
        seq1, ksize), kmer_indexing_brut(seq2, ksize)
    counts_1, counts_2 = Counter({k: v/len(seq1) for k, v in counts_1.items()}), Counter({
        k: v/len(seq2) for k, v in counts_2.items()})
    counts_1.subtract(counts_2)
    figure = plt.figure(figsize=(7, 20))
    plt.plot(counts_1.values(), counts_1.keys())
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.savefig(f"{OUTPUT_PATH}delta_kmers.png", bbox_inches='tight')


if __name__ == "__main__":
    OUTPUT_PATH: str = f"output/figures/compdiff/"
    listing = [f"{a}{b}" for a in ['A', 'T', 'G', 'C']
               for b in ['A', 'T', 'G', 'C']]
    for level in ['domain', 'phylum', 'group', 'order', 'family']:
        elts = compute_signatures(
            level, 'genomes/143_prokaryote_genomes', listing)
        for key, elt in elts.items():
            print("\n"+key+"\n")
            print(elt)
            fig = plt.figure()
            cm = plt.get_cmap('rainbow', 3)
            ax = fig.add_subplot(111)
            cax = ax.matshow(elt, cmap=cm, vmin=0, vmax=2)
            plt.title(f"{key}")
            plt.xticks(rotation=90)
            plt.yticks(fontsize=9)
            plt.xticks(fontsize=9)
            ax.set_xticks([i for i in range(16)])
            ax.set_yticks([i for i in range(16)])
            ax.set_xticklabels([listing[i] for i in range(16)])
            ax.set_yticklabels([listing[i] for i in range(16)])
            ax.tick_params(axis=u'both', which=u'both', length=0)
            #cbar = fig.colorbar(cax)
            #cbar.ax.set_yticks([0, 1, 2])
            #cbar.ax.set_yticklabels(['$f < \mu - \sigma$', '$f = \mu \pm \sigma$', '$f > \mu + \sigma$'])
            plt.savefig(f"{OUTPUT_PATH}{level}/{key}_compdiff_nocbar.png",
                        bbox_inches='tight')

    """
    parser = ArgumentParser()
    # declaring args
    parser.add_argument(
        "name", help="db name", type=str)
    # executing args
    args = parser.parse_args()
    # plot_repartition_top_kmers(6, my_parser("genomes/sequence.fna", True, True, "merge")['merge'], "1111", 4)
    # delta_sequence(my_parser("genomes/sequence.fna", True, True, "merge")['merge'], my_parser("genomes/sequence_2.fna", True, True, "merge")['merge'], "1111", 4)
    Path(OUTPUT_PATH).mkdir(parents=True, exist_ok=True)
    plot_database_features(args.name)
    """
