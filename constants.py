# this file is meant to fetch default values to the software, mostly for support functions

# for utilities.py
from typing import Callable


EMAIL: str = 'siegfried.dubois@inria.fr'
ANNOTATE_PATH: str = "genomes/to_annotate"
SUMMARY_FILE: str = "genomes/assembly_summary.txt"
UNK_PATH: str = "genomes/unk_small/"
TRAIN_PATH: str = "genomes/train_small"

# for wisp.py
DATABASE: str = "tmp_database"
PARAMS: str = "wisp_params_pipeline.json"
SAMPLE_PATH: str = "genomes/unk_small/"
PREFIX_JOB: str = "onevsall"
LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family', 'merged']

# for plotters.py
OUTPUT_PATH: str = "output/figures/"
