# this file is meant to fetch default values to the software, mostly for support functions

# for utilities.py
from typing import Callable
from wisp_lib import kmer_indexing, read_and_its_compl


EMAIL: str = 'siegfried.dubois@inria.fr'
ANNOTATE_PATH: str = "/udd/sidubois/Stage/Genomes/to_annotate"
SUMMARY_FILE: str = "/udd/sidubois/Stage/assembly_summary.txt"
UNK_PATH: str = "genomes/unk/"

# for wisp.py
DATABASE: str = "4k_400genomes_1111"
PARAMS: str = "wisp_params.json"
SAMPLE_PATH: str = "genomes/unk/"
PREFIX_JOB: str = "tex"

# for sample_class.py
COUNT_METHOD: Callable = kmer_indexing
