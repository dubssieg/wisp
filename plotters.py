"Analysis upon kmers repartition and stuff"

import matplotlib.pyplot as plt
import pandas as pd
from wisp_lib import kmer_indexing_brut
from python_tools import my_parser
from constants import OUTPUT_PATH


def plot_repartition_top_kmers(number_to_plot: int, sequence: str, pattern: str, ksize: int) -> None:
    # gives most common at global scale
    counter = kmer_indexing_brut(
        sequence, ksize, pattern).most_common(number_to_plot)
    elements = [e for (e, _) in counter]
    df = pd.DataFrame(columns=elements)
    # for those, we compute locally their abundance
    for i in range(0, len(sequence), int(len(sequence)/50000)):
        local_kmers = kmer_indexing_brut(sequence[i:i+50000], ksize, pattern)
        my_local_kmers = {k: v for k,
                          v in local_kmers.items() if k in elements}
        label = f"[{i}:{i+50000}]"
        my_series = pd.Series(data=my_local_kmers, index=list(
            my_local_kmers.keys()), name=label)
        df = df.append(my_series).fillna(0)
    df = df.transpose()
    ax = df.plot(figsize=(
        20, 6), ylabel=f'kmers/{int(len(sequence)/5000)}bp', rot=90, colormap='cividis')
    plt.savefig(f"{OUTPUT_PATH}kmers_repartition", bbox_inches='tight')


def delta_sequence(seq1: str, seq2: str, pattern: str, ksize: int) -> None:
    counts_1, counts_2 = kmer_indexing_brut(
        seq1, ksize, pattern), kmer_indexing_brut(seq2, ksize, pattern)
    plt.plot(counts_1.subtract(counts_2))
    plt.savefig(f"{OUTPUT_PATH}delta_kmers", bbox_inches='tight')


#plot_repartition_top_kmers(6, my_parser("genomes/sequence.fna", True, True, "merge")['merge'], "1111", 4)

delta_sequence(my_parser("genomes/sequence.fna", True, True, "merge")
               ['merge'], my_parser("genomes/sequence_2.fna", True, True, "merge")['merge'], "1111", 4)
