from cmath import nan
from python_tools import my_parser
import matplotlib.pyplot as plt
from wisp_lib import decode_kmer_4
from statistics import mean
from wisp_lib import species_list
from sample_class import make_datasets
from build import make_model
from main import test_model


def bacterial_signature(sp: str, func=mean, ratio=3/2) -> dict:
    """
    Create a bacterial signature of the given genome based upon most frequent kmers

    * seq_path : complete path to genome
    * func : selection method for threshold
    * ratio : multiplier for func(values)
    """
    SEQUENCE = f"/udd/sidubois/Stage/Genomes/{sp}.fna"
    parsed_gen = my_parser(SEQUENCE, True, True, f"{sp}")
    indexed_gen = indexing_by_signature_with_subsampling_linear(parsed_gen)
    #values = [float(v) for v in indexed_gen[f"{sp}"].values()]
    # if float(v) > func(values)*(ratio)
    print(indexed_gen)
    return {decode_kmer_4(k): float(v) for k, v in indexed_gen[f"{sp}"].items()}


def plotter(list_signatures: list) -> None:
    """
    Create a representation for a list of signatures, replacing absent valuses by NaN.

    * list_signatures : a list of str+dicts, from bacterial_signature
    """
    set_keys: list = []
    for spc in list_signatures:
        set_keys = [*set_keys, *spc[1].keys()]
    sorted(set(set_keys))

    for dico in list_signatures:
        y_plot = [dico[1][k] if k in dico[1].keys() else float('NaN')
                  for k in set_keys]
        plt.scatter(set_keys, y_plot, label=dico[0], marker='o')

    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc='lower left',
               mode='expand', ncol=len(list_signatures))
    plt.show()


def plotter_bars_grouped() -> None:
    """
    Create a representation for a list of signatures, replacing absent valuses by NaN.

    * list_signatures : a list of str+dicts, from bacterial_signature
    """
    INPUT_PATH = '/udd/sidubois/Stage/Genomes/'
    OUTPUT_PATH = 'data/'
    JOB = 'temp'

    graph_values: list = []
    graph_keys: list = []

    for km in [4, 5, 6]:
        for fnc in [(None, 0), (mean, 1.2), (mean, 1.3), (mean, 1.4), (mean, 1.5)]:

            make_datasets(INPUT_PATH, OUTPUT_PATH, JOB,
                          func=fnc[0], ratio=fnc[1], kmer_size=km)
            make_model(OUTPUT_PATH, JOB)
            graph_values.append(test_model(OUTPUT_PATH, JOB))
            graph_keys.append(
                f"{km}mer {'w/o.' if fnc[0]==None else 'w. '+str(fnc[1])}")

    y_pos = range(len(graph_keys))
    plt.bar(graph_keys, graph_values)
    plt.xticks(y_pos, graph_keys, rotation=90)
    plt.show()


if __name__ == "__main__":
    "Executes main procedure"
    plotter_bars_grouped()
    #species = species_list('/udd/sidubois/Stage/Genomes/')
    #plotter([(sp, bacterial_signature(sp)) for sp in species])

# IDEE : peut on extraire une "signature" de la bactérie en prenant tout ce qui est au delà d'un seuil? 3eme quartile + ?
