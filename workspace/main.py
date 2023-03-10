#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from os import walk, path
from rich.traceback import install
from rich import print
from tharospytools import futures_collector
from workspace.create_database import build_database
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
    help="Input folder containig genomes",
    type=str
)

## Subparser for model creation ##

parser_model: ArgumentParser = subparsers.add_parser(
    'model',
    help="Creates the model from the specified database.")

parser_model.add_argument(
    "database_path",
    help="Path for database",
    type=str
)

parser_model.add_argument(
    "-p",
    "--parameters",
    help="Specifies a parameter file",
    type=str,
    default=f'{path.dirname(__file__)}/parameters_files/params.json'
)

#######################################

args = parser.parse_args()


def main() -> None:
    "Main call for subprograms"
    if len(argv) == 1:
        print(
            "[dark_orange]You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands."
        )
        exit()

    if args.subcommands == 'build':
        print("[dark_orange]Starting database creation")
        output_path: str = build_database(args.parameters, args.database_name,
                                          [path.abspath(path.join(dirpath, f)) for dirpath, _, filenames in walk(args.input_folder) for f in filenames])
        print(
            f"[dark_orange]Database sucessfully built @ {output_path}")

    if args.subcommands == 'model':
        pass
