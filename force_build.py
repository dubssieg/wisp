from sample_class import make_datasets
from build_softprob import make_model, init_parameters
from os import listdir
from python_tools import my_function_timer, my_output_msg, my_logs_global_config, my_logs_clear
from wisp_lib import load_mapping, load_json
from argparse import ArgumentParser
from constants import FUNC, RATIO, TAXAS_LEVELS


@my_function_timer("Building full database")
def build_full_db(args) -> None:
    """For export or pre-processing purposes

    Args:
        db_name (str): _description_

    Raises:
        ValueError: _description_
    """

    # we try to load parmas file and gather data from it
    try:
        my_params: dict = load_json(args.params)
        # storing args
        JOB: str = "build_full_db"
        DATABASE: str = args.database_name
        INPUT_PATH: str = my_params['input']
        OUTPUT_PATH: str = my_params['output']
        nr = int(my_params['nb_boosts'])
        tree_depth = int(my_params['tree_depth'])
    # if any error happens
    except:
        raise ValueError(
            "Incorrect or missing parameters file ; check path and/or contents of json reference.")

    list_of_genomes = [genome.split('.')[0]
                       for genome in listdir(f"{INPUT_PATH}train/")]
    for taxa in TAXAS_LEVELS:
        KMER_SIZE_REF, RS_REF, SAMPLING_REF, PATTERN_REF = my_params[f"{taxa}_ref"]

        list_parent_level = [i for i in set([e.split('_')[TAXAS_LEVELS.index(
            taxa)-1] for e in list_of_genomes])] if taxa != 'domain' else [False]

        for parent_level in list_parent_level:
            if isinstance(parent_level, bool):
                parent_level = None

            my_output_msg(
                f"Building dataset at level {taxa} for parent level {parent_level}")

            make_datasets(
                input_style=False,
                job_name=JOB,
                input_dir=INPUT_PATH,
                path=OUTPUT_PATH,
                datas=['train', 'test'],
                db_name=DATABASE,
                sampling=SAMPLING_REF,
                kmer_size=KMER_SIZE_REF,
                func=FUNC,
                ratio=RATIO,
                read_size=RS_REF,
                classif_level=taxa,
                sp_determied=parent_level,
                pattern=PATTERN_REF
            )

            map_sp = load_mapping(OUTPUT_PATH, DATABASE,
                                  taxa, parent_level)

            make_model(JOB, OUTPUT_PATH, taxa, DATABASE,
                       parent_level, init_parameters(len(map_sp), tree_depth), number_rounds=nr)


if __name__ == "__main__":
    "Executes main procedure"
    my_logs_global_config("LOG_wisp")
    my_logs_clear("LOG_wisp.log")  # we clean before entering loop
    parser = ArgumentParser()

    # declaring args
    parser.add_argument(
        "database_name", help="name of database, builded if not exists", type=str)
    parser.add_argument(
        "params", help="path to a params .json file", type=str)

    # executing args
    args = parser.parse_args()

    build_full_db(args)
