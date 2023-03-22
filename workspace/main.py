#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from os import walk, path
from treelib import Tree
from rich.traceback import install
from rich import print
from tharospytools import futures_collector
from json import load
from workspace.create_database import build_database
from workspace.create_model import make_model
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
        exit()

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

        print(phylo_tree)

        nodes_per_level: dict = {level: [node.tag for node in list(phylo_tree.filter_nodes(
            lambda x: phylo_tree.depth(x) == i))] for i, level in enumerate(['domain', 'phylum', 'group', 'order'], start=1)}

        print(
            "[dark_orange]Starting model creation"
        )

        # Loading data => should be put in the main call to escape loading it at each iteration
        with open(output_path, 'r', encoding='utf-8') as jdb:
            datas = load(jdb)

        retcodes: list = futures_collector(make_model, [(datas, output_path, taxonomic_level, target_taxa)
                                                        for taxonomic_level, targets in nodes_per_level.items() for target_taxa in targets])

        for retcode in retcodes:
            print(
                f"[dark_orange]Model sucessfully built @ {retcode}"
            )
