from argparse import ArgumentParser
from sample_class import make_datasets
from build_softprob import make_model, init_parameters, make_testing
from warnings import filterwarnings
from python_tools import my_logs_global_config
from datetime import datetime
from wisp_view import gen_html_report, tree_render, plot_boosting
from wisp_lib import check_if_database_exists, check_if_model_exists, load_mapping, load_json
from predictors import test_unk_sample, save_output, test_model
from constants import RATIO, FUNC, TAXAS_LEVELS
from pathlib import Path

if __name__ == "__main__":
    "Executes main procedure"
    filterwarnings('ignore')  # to ignore xgboost warnnings
    my_logs_global_config("LOG_wisp")

    parser = ArgumentParser()

    # declaring args
    parser.add_argument(
        "database_name", help="name of database, builded if not exists", type=str)
    parser.add_argument(
        "params", help="path to a params .json file", type=str)
    parser.add_argument("job_name", help="name of the given job", type=str)
    parser.add_argument(
        "-f", "--file", type=str, help="file to load. if not specified, loads the full unk/ folder")

    # executing args
    args = parser.parse_args()

    input_file: str | bool = args.file if args.file != None else False

    # we try to load parmas file and gather data from it
    try:
        my_params: dict = load_json(args.params)
        # storing args
        JOB: str = args.job_name
        DATABASE: str = args.database_name
        INPUT_PATH: str = my_params['input']
        OUTPUT_PATH: str = my_params['output']
        nr = int(my_params['nb_boosts'])
        threshold = float(my_params['threshold'])
        reads_threshold = float(my_params['reads_th'])
        test_state = bool(my_params['full_test_set'])
        force_rebuild = bool(my_params['force_model_rebuild'])
        tree_depth = int(my_params['tree_depth'])
        func_reads = str(my_params['selection_mode'])
    # if any error happens
    except:
        raise ValueError(
            "Incorrect or missing parameters file ; check path and/or contents of json reference.")

    # init iterables and memory spaces
    topmost: dict = {}
    test_results: dict = {}
    output: dict = {'domain': None, 'phylum': None,
                    'group': None, 'order': None, 'family': None}
    Path(f"output/{JOB}/").mkdir(parents=True, exist_ok=True)
    for taxa in TAXAS_LEVELS:
        KMER_SIZE_REF, RS_REF, SAMPLING_REF = my_params[f"{taxa}_ref"]
        KMER_SIZE_SAMPLE, RS_SAMPLE, SAMPLING_SAMPLE = my_params[f"{taxa}_sample"]

        list_parent_level = output[
            f"Possible for {TAXAS_LEVELS[TAXAS_LEVELS.index(taxa)-1]}"] if taxa != 'domain' else [False]

        for parent_level in list_parent_level:
            if isinstance(parent_level, bool):
                parent_level = None

            if not check_if_database_exists(DATABASE, OUTPUT_PATH, taxa, parent_level):

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
                    sp_determied=parent_level
                )

            map_sp = load_mapping(OUTPUT_PATH, DATABASE,
                                  taxa, parent_level)

            if force_rebuild or not check_if_model_exists(DATABASE, OUTPUT_PATH, taxa, parent_level):

                make_model(JOB, OUTPUT_PATH, taxa, DATABASE,
                           parent_level, init_parameters(len(map_sp), tree_depth), number_rounds=nr)

            make_datasets(
                input_style=input_file,
                job_name=JOB,
                input_dir=INPUT_PATH,
                path=OUTPUT_PATH,
                datas=['unk'],
                db_name=DATABASE,
                sampling=SAMPLING_SAMPLE,
                kmer_size=KMER_SIZE_SAMPLE,
                func=FUNC,
                ratio=RATIO,
                read_size=RS_SAMPLE,
                classif_level=taxa,
                sp_determied=parent_level
            )

            # full test set, takes time, but gives info on structure
            if test_state:
                successive_boost_results = make_testing(
                    size_kmer=KMER_SIZE_REF,
                    job_name=JOB,
                    sp_determined=parent_level,
                    path=OUTPUT_PATH,
                    db_name=DATABASE,
                    classif_level=taxa,
                    class_count=len(map_sp),
                    model_parameters=init_parameters(len(map_sp)),
                    number_rounds=nr
                )

                plot_boosting(successive_boost_results,
                              JOB, taxa, parent_level, nr)

            # base tests for heatmap and evaluators
            test_results[f"{taxa}_{parent_level}"] = (test_model(
                OUTPUT_PATH, JOB, DATABASE, taxa, reads_threshold, parent_level))

            output_temp = test_unk_sample(
                OUTPUT_PATH, JOB, DATABASE, taxa, parent_level, threshold, reads_threshold, test_state)
            topmost[f"{taxa}_{parent_level}"] = output_temp[f"Reads summation {taxa}"]

            if f"Possible for {taxa}" in output:
                output[f"Possible for {taxa}"] = [
                    *output[f"Possible for {taxa}"], *output_temp[f"Possible for {taxa}"]]

            output = {**output_temp, **output}

    # extraction of most probable path
    output["domain"] = topmost[f"domain_None"]
    output["phylum"] = topmost[f"phylum_{output['domain']}"]
    output["group"] = topmost[f"group_{output['phylum']}"]
    output["order"] = topmost[f"order_{output['group']}"]
    output["family"] = topmost[f"family_{output['order']}"]

    path_taxa = [f"{output['domain']} (d)", f"{output['phylum']} (p)",
                 f"{output['group']} (g)", f"{output['order']} (o)", f"{output['family']} (f)"]

    output["parcimonious_path"] = tree_render(output, JOB, path_taxa)

    save_output({'Date': f"{datetime.today().strftime('%Y.%m.%d - %H:%M:%S')}", **
                 vars(args), **output}, JOB)
    gen_html_report(my_params, JOB, [], output, TAXAS_LEVELS,
                    test_results, threshold, test_state, round(reads_threshold, 2))
