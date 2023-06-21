#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from os import walk, path
from json import load, dump
from pickle import dump as pdump, load as pload
from pathlib import Path
from rich.traceback import install
from rich import print
from treelib import Tree
from Bio import SeqIO
from tharospytools import futures_collector
from workspace.create_database import build_database
from workspace.create_model import make_model
from workspace.create_prediction import prediction


parser: ArgumentParser = ArgumentParser(
    description='Bacteria family identification tool.', add_help=True)
parser.add_argument(
    "-l", "--locals", help="Display locals on error.", action='store_true')
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

parser_prediction: ArgumentParser = subparsers.add_parser(
    'predict',
    help="Creates the samples and evaluates them.")

parser_prediction.add_argument(
    "-p",
    "--parameters",
    help="Specifies a parameter file",
    type=str,
    default=f'{path.dirname(__file__)}/parameters_files/params.json'
)

parser_prediction.add_argument(
    "database_name",
    help="Name for database",
    type=str
)

parser_prediction.add_argument(
    "input_folder",
    help="Input folder containig unknown genomes",
    type=str
)

parser_prediction.add_argument(
    "output_folder",
    help="Input folder containig reference genomes",
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

    install(show_locals=args.locals)

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

        phylo_tree.show()  # data_property='code'

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
            _, _, taxonomic_level, target_taxa = fargs[i]

            if model_path is not None and config_path is not None:

                try:
                    node = phylo_tree[f"{target_taxa.lower()}_{taxonomic_level}"]
                    node.data.model_path = model_path
                    node.data.config_path = config_path
                except KeyError:
                    phylo_tree.remove_node(target_taxa.lower())

        with open(phylo_path := f"{path.dirname(__file__)}/model/{args.database_name}_phylo_tree.txt", 'wb') as jtree:
            pdump(phylo_tree, jtree)

        print(f"[dark_orange]Finished computing models, tree @ {phylo_path}")
        exit(0)

    if args.subcommands == 'predict':
        print(
            "[dark_orange]Starting sample creation"
        )
        # Loading params file
        with open(args.parameters, 'r', encoding='utf-8') as pfile:
            params: dict = load(pfile)

        # Loading the phylogenetic tree
        with open(phylo_path := f"{path.dirname(__file__)}/model/{args.database_name}_phylo_tree.txt", 'rb') as jtree:
            tree: Tree = pload(jtree)

        try:
            threshold: float = params["threshold"]
            if threshold > 1.0:
                threshold = 1.0
            elif threshold < 0.01:
                threshold = 0.01
        except KeyError as exc:
            raise RuntimeError(
                "Invalid parameter file, must contain a read acceptance threshold value between 0.01 (1% identity) and 1.0 (100% identity).") from exc

        # iterating over input files
        for genome in [path.abspath(path.join(dirpath, f)) for dirpath, _, filenames in walk(
                args.input_folder) for f in filenames]:

            with open(genome, 'r', encoding='utf-8') as freader:
                genome_data: dict = {fasta.id: str(fasta.seq)
                                     for fasta in SeqIO.parse(freader, 'fasta')}

            prediction_results: list = futures_collector(prediction, pargs := [
                (id_sequence, dna_sequence, params, tree, threshold, args.parameters) for id_sequence, dna_sequence in genome_data.items()])

            with open(report_path := path.join(args.output_folder, f"{Path(genome).stem}_job_output.json"), 'w', encoding='utf-8') as jwriter:
                dump({identifier: prediction_results[i] for i, (identifier, _, _, _, _, _) in enumerate(
                    pargs)}, jwriter)

            print(
                f"[dark_orange]Job on file {Path(genome).stem} ended sucessfully, report @ {report_path}"
            )
        exit(0)
