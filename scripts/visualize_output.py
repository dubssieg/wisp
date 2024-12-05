"""Aims to display results in a Sankey graph"""
from pickle import load as pload
from json import load as jload
from argparse import ArgumentParser, SUPPRESS
from pathlib import Path
from os import path
from random import choice
from rich.traceback import install
from treelib import Tree
from tharospytools.matplotlib_tools import get_palette
from plotly import graph_objects as go
from dash import Dash, dcc, html
from warnings import filterwarnings


def dash_app(fig, fig2, job_name: str = 'Job report'):

    app: Dash = Dash()
    app.layout = html.Div(
        [html.H1(
            children=job_name,
            style={
                'textAlign': 'center',
                'color': 'deepslateblue'
            }
        ),
            # html.Div([dcc.Graph(figure=fig)], style={'display': 'inline-block', 'width': '100%'}),
            html.Div([
                html.Div([
                    dcc.Graph(figure=fig)], style={'display': 'inline-block', 'width': '70%'}
                ),
                html.Div([
                    dcc.Graph(figure=fig2)], style={'display': 'inline-block', 'width': '30%'}
                )
            ], style={'align': 'center'})
        ]
    )

    app.run_server(debug=True, use_reloader=True)


def plot_report(phylo_path: str, json_report: str, output_path: str, reads_per_sample: int = 100, showing_rejected: bool = False, cutoff: int = 50) -> None:
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
        len(output_data)*2, cmap_name='Greens', as_hex=True)[:: -1]
    for idx_col, col in enumerate(palette_links):
        r, g, b = tuple(int(col[1:][i: i+2], 16) for i in (0, 2, 4))
        palette_links[idx_col] = f"rgba({r},{g},{b},0.2)"

    for idx_sample, (sequence, values) in enumerate(output_data.items()):
        if values == "REJECTED" and showing_rejected:
            seqnames.append(sequence)
            source.append(0)
            target.append(len(labels)-1)
            value.append(reads_per_sample)
            color_links.append('black')
        elif values != 'REJECTED':
            for idx_level, level in enumerate(values):
                for taxa_level, results_per_level in level.items():
                    for taxunum, count in results_per_level.items():
                        if count >= cutoff:
                            if f"{taxunum}_{idx_level+1}" in mappings_taxa:
                                seqnames.append(sequence)
                                source.append(
                                    mappings_taxa[f"{taxa_level}_{idx_level}"])
                                target.append(
                                    mappings_taxa[f"{taxunum}_{idx_level+1}"])
                                value.append(count)
                                color_links.append(choice(palette_links))

    sunburst_counts: dict = dict()
    for idx, v in enumerate(value):
        if target[idx] in inv_map:
            if inv_map[target[idx]] in sunburst_counts:
                sunburst_counts[inv_map[target[idx]]] += v
            else:
                sunburst_counts[inv_map[target[idx]]] = v
        else:
            if 'REJECTED' in sunburst_counts:
                sunburst_counts['REJECTED'] += v
            else:
                sunburst_counts['REJECTED'] = v

    print(sunburst_counts, file=open(path.join(output_path, f"{Path(json_report).stem}_counts.json"), 'w', encoding='utf-8'))

    value = [(v/reads_per_sample)*100 for v in value]

    data_sankey = go.Sankey(link=dict(source=source, target=target, value=value, label=seqnames, color=color_links), node=dict(label=labels, color=color), valueformat=".0f",
                            valuesuffix="%")

    mappings_ancestors: dict = {
        f"{node.tag}_{tree.depth(node)}": f"{tree.get_node(tree.ancestor(node.identifier)).tag}_{tree.depth(node)-1}" for i, node in enumerate(tree.all_nodes_itr()) if tree.depth(node) > 0
    }
    mappings_ancestors["REJECTED"] = "Root_0"
    sunburst_points: list[str] = [
        x for x in mappings_ancestors.keys()]
    sunburst_parents: list[str] = [mappings_ancestors[x]
                                   for x in sunburst_points]
    sunburst_values: list[int] = [sunburst_counts[x]
                                  if x in sunburst_counts else 0 for x in sunburst_points]
    sunburst_labels: list[str] = [
        f"{label.split('_')[0]}<br>{sunburst_values[i]}" for i, label in enumerate(sunburst_points)]

    data_sunburst = go.Sunburst(ids=sunburst_points, labels=sunburst_labels,
                                parents=sunburst_parents, values=sunburst_values, insidetextorientation='radial')

    layout = go.Layout(
        autosize=False,
        width=800,
        height=2000
    )
    layout2 = go.Layout(
        autosize=False,
        width=1000,
        height=1000
    )

    fig = go.Figure(data_sankey, layout=layout)

    fig2 = go.Figure(data_sunburst, layout=layout2)

    fig.update_layout(
        margin=dict(t=0, l=0, r=0, b=0),
        font=dict(
            # family="Courier New, monospace",
            size=8,  # Set the font size here
            color="RebeccaPurple"
        )
    )
    fig2.update_layout(
        margin=dict(t=0, l=0, r=0, b=0),
        font=dict(
            # family="Courier New, monospace",
            size=8,  # Set the font size here
            color="RebeccaPurple"
        )
    )

    fig.write_image(
        path.join(output_path, f"{Path(json_report).stem}_sankey_diag.svg"))
    fig2.write_image(
        path.join(output_path, f"{Path(json_report).stem}_sunburst_diag.svg"))

    # dash_app(fig, fig2)

    # fig.write_html(path.join(output_path, f"{Path(json_report).stem}_graph.html"))


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
    filterwarnings("ignore", category=DeprecationWarning)

    plot_report(args.taxonomy, args.report, args.output)
