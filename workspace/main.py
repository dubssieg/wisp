#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from os import walk, path
from json import load
from rich.traceback import install
from rich import print
from tharospytools import futures_collector
from workspace.create_database import build_database
from workspace.create_model import make_model
from workspace.create_sample import build_sample
from workspace.create_prediction import make_prediction
install(show_locals=True)

parser: ArgumentParser = ArgumentParser(
    description='Bacteria family identification tool.', add_help=True)
subparsers = parser.add_subparsers(
    help='Available subcommands', dest="subcommands")

parser._positionals.title = 'Subcommands'
parser._optionals.title = 'Global Arguments'

## Subparser for database creation ##

parser_database: ArgumentParser = subparsers.add_parser(
    'build',
    help="Creates the database from the specified set of files.")

parser_database.add_argument(
    "-p",
    "--parameters",
    help="Specifies a parameter file",
    type=str,
    default=f'{path.dirname(__file__)}/parameters_files/params.json'
)

parser_database.add_argument(
    "database_name",
    help="Name for database",
    type=str
)

parser_database.add_argument(
    "input_folder",
    help="Input folder containig reference genomes",
    type=str
)

#######################################

## Subparser for prediction ##

parser_database: ArgumentParser = subparsers.add_parser(
    'predict',
    help="Creates the samples and evaluates them.")

parser_database.add_argument(
    "-p",
    "--parameters",
    help="Specifies a parameter file",
    type=str,
    default=f'{path.dirname(__file__)}/parameters_files/params.json'
)

parser_database.add_argument(
    "database_name",
    help="Name for database",
    type=str
)

parser_database.add_argument(
    "input_folder",
    help="Input folder containig unknown genomes",
    type=str
)

#######################################

args = parser.parse_args()


def main() -> None:
    "Main call for subprograms"
    if len(argv) == 1:
        print(
            "[red]You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands."
        )
        exit(2)

    if args.subcommands == 'build':
        print(
            "[dark_orange]Starting database creation"
        )
        output_path, phylo_tree = build_database(
            args.parameters,
            args.database_name,
            [path.abspath(path.join(dirpath, f)) for dirpath, _, filenames in walk(
                args.input_folder) for f in filenames]
        )
        print(
            f"[dark_orange]Database sucessfully built @ {output_path}"
        )

        phylo_tree.show()

        nodes_per_level: dict = {level: [node.tag for node in list(phylo_tree.filter_nodes(
            lambda x: phylo_tree.depth(x) == i))] for i, level in enumerate(['root', 'domain', 'phylum', 'group', 'order'], start=0)}

        print(
            "[dark_orange]Starting model creation"
        )

        # Loading data => should be put in the main call to escape loading it at each iteration
        with open(output_path, 'r', encoding='utf-8') as jdb:
            datas = load(jdb)

        retcodes: list = futures_collector(make_model, fargs := [(datas, output_path, taxonomic_level, target_taxa)
                                                                 for taxonomic_level, targets in nodes_per_level.items() for target_taxa in targets])

        for i, (model_path, config_path) in enumerate(retcodes):
            _, _, _, target_taxa = fargs[i]

            try:
                node = phylo_tree[target_taxa.lower()]
                node.data.model_path = model_path
                node.data.config_path = config_path
            except KeyError:
                pass

            # print(f"[dark_orange]Model for {target_taxa} (level {taxonomic_level}) sucessfully built @ {model_path}")

        phylo_tree.show(data_property="model_path")
        exit(0)

    if args.subcommands == 'predict':
        print(
            "[red]Currently non fully implemented."
        )
        print(
            "[dark_orange]Starting database creation"
        )
        output_path = build_sample(
            args.parameters,
            [path.abspath(path.join(dirpath, f)) for dirpath, _, filenames in walk(
                args.input_folder) for f in filenames]
        )
        print(
            f"[dark_orange]Query sucessfully built @ {output_path}"
        )
        exit(126)
        results = make_prediction(
            "/usr/local/lib/python3.11/site-packages/workspace/model/test/Lactobacillales.json",
            "/usr/local/lib/python3.11/site-packages/workspace/model/test/Lactobacillales_params.json",
            output_path,
            normalisation_func='delta_mean',
            read_identity_threshold=0.8
        )
        print(results)

        # Loading the phylogenetic tree

        # Loading the model in a booster

        # Evaluate at one level
        # Use the tree to select next level
