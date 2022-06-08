# this file is meant to fetch default values to the software, mostly for support functions

# for utilities.py
from typing import Callable


EMAIL: str = 'siegfried.dubois@inria.fr'
ANNOTATE_PATH: str = "genomes/to_annotate"
SUMMARY_FILE: str = "genomes/assembly_summary.txt"
UNK_PATH: str = "genomes/unk/"

# for wisp.py
DATABASE: str = "129genomes"
PARAMS: str = "wisp_params.json"
SAMPLE_PATH: str = "genomes/unk/"
PREFIX_JOB: str = "tex"
LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family', 'merged']

# for plotters.py
OUTPUT_PATH: str = "output/figures/"
