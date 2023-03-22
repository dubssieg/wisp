from sample_class import make_datasets
from build_softprob import make_model, init_parameters
from os import listdir
from wisp_tools import my_function_timer, my_output_msg, my_logs_global_config, my_logs_clear
from wisp_lib import load_mapping, load_json, check_if_database_exists, check_if_merged_model_exists, check_if_model_exists, check_if_merged_database_exists

from traceback import format_exc


@my_function_timer("Building full database")
def build_full_db(args) -> None:
    """Builds full database, step by step, with all its models.

    Args:
        args (Namespace): args of call

    Raises:
        ValueError: _description_
    """
    # we try to load params file and gather data from it
    try:
        my_params: dict = load_json(args.params)
        # storing args
        JOB: str = "build_full_db"
        DATABASE: str = args.database_name
        TRAIN_PATH: str = my_params['input_train']
        DATABASE_PATH: str = my_params['database_output']
        TAXAS_LEVELS: list[str] = [
            t for t in my_params['levels_list'] if t in args.levels]
        GLOBAL_TAXAS_LEVELS: list[str] = my_params['levels_list']
        nr: int = int(my_params['nb_boosts'])
        tree_depth: int = int(my_params['tree_depth'])
        force_rebuild: bool = bool(my_params['force_model_rebuild'])
        WINDOW: int = my_params['window_size']
        KMER_SIZE_MERGED_REF, SAMPLING_MERGED_REF, PATTERN_MERGED_REF = my_params[
            f"merged_ref"]
    # if any error happens
    except Exception as exc:
        my_output_msg(format_exc())
        raise ValueError(
            "Incorrect or missing parameters file ; check path and/or contents of json reference.") from exc

    list_of_genomes = [genome.split('.')[0]
                       for genome in listdir(f"{TRAIN_PATH}")]
    for taxa in TAXAS_LEVELS:
        if taxa != 'merged':
            KMER_SIZE_REF, SAMPLING_REF, PATTERN_REF = my_params[f"{taxa}_ref"]

            list_parent_level = [i for i in set([e.split('_')[GLOBAL_TAXAS_LEVELS.index(
                taxa)-1] for e in list_of_genomes])] if taxa != 'domain' else [False]

            for parent_level in list_parent_level:
                if isinstance(parent_level, bool):
                    parent_level = None

                if not check_if_database_exists(DATABASE, DATABASE_PATH, taxa, parent_level):

                    make_datasets(
                        job_name=JOB,
                        input_dir=TRAIN_PATH,
                        path=DATABASE_PATH,
                        datas=['train', 'test'],
                        db_name=DATABASE,
                        sampling=SAMPLING_REF,
                        kmer_size=KMER_SIZE_REF,
                        read_size=WINDOW,
                        classif_level=taxa,
                        sp_determied=parent_level,
                        pattern=PATTERN_REF
                    )
