# this file is meant to fetch default values to the software, mostly for support functions

# for utilities.py
EMAIL: str = 'siegfried.dubois@inria.fr'
ANNOTATE_PATH: str = "/udd/sidubois/Stage/143_Genomes/train"
SUMMARY_FILE: str = "/udd/sidubois/Stage/assembly_summary.txt"

# old parameters for main.py
FUNC: str | None = None
RATIO: float = 1.5

# iterables
TAXAS_LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']

# for wisp.py
DATABASE: str = "4k_143genomes_1111"
PARAMS: str = "wisp_params.json"
SAMPLE_PATH: str = "/udd/sidubois/Stage/Genomes/unk/"
PREFIX_JOB: str = "1111_143genomes_deltamean_025"
THREADPOOL: int = 8

# for mass_analysis.py
FOLDER_LIST: list[str] = ['purge_reads_0,00_8_boosts',
                          'purge_reads_0,00_10_boosts',
                          'purge_reads_0,50_8_boosts',
                          'purge_reads_0,50_10_boosts',
                          '2_purge_reads_0,00_10_boosts',
                          '2_purge_reads_0,25_10_boosts',
                          'purge_reads_0,25_8_boosts',
                          'purge_reads_0,25_10_boosts',
                          'purge_reads_0,25_12_boosts']
