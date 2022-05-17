# this file is meant to fetch default values to the software, mostly for support functions

# for script utilities.py
ANNOTATE_PATH: str = "/udd/sidubois/Stage/Genomes/to_annotate"
SUMMARY_FILE: str = "/udd/sidubois/Stage/assembly_summary.txt"

# old parameters for main
FUNC: str | None = None
RATIO: float = 1.5

# iterables
TAXAS_LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']

# for wisp.py
DATABASE: str = "4k_16_05_2022"
PARAMS: str = "wisp_params.json"
SAMPLE_PATH: str = "/udd/sidubois/Stage/Genomes/unk/"
PREFIX_JOB: str = "mda"
