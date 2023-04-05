"Aims to display results in a Sankey graph"
from pickle import load as pload
from json import load as jload
from argparse import ArgumentParser, SUPPRESS
from rich.traceback import install
from treelib import Tree
from pathlib import Path
from os import path
from tharospytools import get_palette
import plotly.graph_objects as go


def plot_report(phylo_path: str, json_report: str, output_path: str, reads_per_sample: int = 500) -> None:
    """Plots a interactive diagram for reads

    Args:
        phylo_path (str): _description_
        json_report (str): _description_
        output_path (str): _description_
        reads_per_sample (int, optional): _description_. Defaults to 500.
    """

    with open(phylo_path, 'rb') as jtree:
        tree: Tree = pload(jtree)

    mappings_taxa: dict = {
        f"{node.tag}_{tree.depth(node)}": i for i, node in enumerate(tree.all_nodes_itr())}
    inv_map = {v: k for k, v in mappings_taxa.items()}

    with open(json_report, 'r', encoding='utf-8') as jreader:
        output_data = jload(jreader)

    source = []  # points sources
    target = []  # points cibles
    value = []  # valeurs de fin de chaîne
    # labels à ajouter sur les segments
    palette_nodes = get_palette(7, as_hex=True)
    labels = [inv_map[x].split('_')[0]
              for x in range(len(mappings_taxa))] + ['REJECTED']
    color = [palette_nodes[int(inv_map[x].split('_')[1])]
             for x in range(len(mappings_taxa))] + ['black']
    seqnames: list = list()
    color_links: list = list()

    palette_links = get_palette(
        len(output_data)*2, cmap_name='Greens', as_hex=True)[::-1]
    for idx_col, col in enumerate(palette_links):
        r, g, b = tuple(int(col[1:][i:i+2], 16) for i in (0, 2, 4))
        palette_links[idx_col] = f"rgba({r},{g},{b},0.2)"

    for idx_sample, (sequence, values) in enumerate(output_data.items()):
        if values == "REJECTED":
            seqnames.append(sequence)
            source.append(0)
            target.append(len(labels)-1)
            value.append(reads_per_sample)
            color_links.append(palette_links[idx_sample])
        else:
            for idx_level, level in enumerate(values):
                for taxa_level, results_per_level in level.items():
                    for taxunum, count in results_per_level.items():
                        if f"{taxunum}_{idx_level+1}" in mappings_taxa:
                            seqnames.append(sequence)
                            source.append(
                                mappings_taxa[f"{taxa_level}_{idx_level}"])
                            target.append(
                                mappings_taxa[f"{taxunum}_{idx_level+1}"])
                            value.append(count)
                            color_links.append(palette_links[idx_sample])

    value = [(v/reads_per_sample)*100 for v in value]

    data = go.Sankey(link=dict(source=source, target=target, value=value, label=seqnames, color=color_links), node=dict(label=labels, color=color), valueformat=".0f",
                     valuesuffix="%")

    fig = go.Figure(data)

    fig.write_html(
        path.join(output_path, f"{Path(json_report).stem}_graph.html"))


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "taxonomy", type=str, help="Path to wisp taxonomy tree.")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Creates a html report of reads repartition')
    parser.add_argument(
        "-o", "--output", help="Specify a output folder", required=True)
    parser.add_argument(
        "-r", "--report", help="Report to be plotted", required=True)
    args = parser.parse_args()

    install(show_locals=True)

    plot_report(args.taxonomy, args.report, args.output)
