from argparse import ArgumentParser
from traceback import format_exc
from sample_class import make_datasets
from build_softprob import make_model, init_parameters, make_testing
from warnings import filterwarnings
from python_tools import my_logs_global_config, my_output_msg, my_fasta_parser
from datetime import datetime
from wisp_view import gen_html_report, tree_render, plot_boosting, plot_pie_merge, make_doc
from wisp_lib import check_if_database_exists, check_if_model_exists, check_if_merged_database_exists, load_mapping, load_json, check_if_merged_model_exists
from predictors import test_unk_sample, save_output, test_model
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

############################################ LOADING STUFF ###############################################

    # parameter for filepath
    input_file: str | bool = args.file if args.file != None else False

    # fasta_reads: set[str] = my_fasta_parser(f"{INPUT_PATH}{input_file}") # TODO file location error

    # we try to load params file and gather data from it
    try:
        my_params: dict = load_json(args.params)
        # storing args
        JOB: str = args.job_name
        DATABASE: str = args.database_name
        TAXAS_LEVELS: list[str] = my_params['levels_list']
        INPUT_PATH: str = my_params['input']
        OUTPUT_PATH: str = my_params['output']
        nr = int(my_params['nb_boosts'])
        threshold = float(my_params['threshold'])
        reads_threshold = float(my_params['reads_th'])
        test_state = bool(my_params['full_test_set'])
        force_rebuild = bool(my_params['force_model_rebuild'])
        tree_depth = int(my_params['tree_depth'])
        func_reads = str(my_params['selection_mode'])
        single_way = bool(my_params['single_way'])
        targeted_level = str(my_params['targeted_level'])
        KMER_SIZE_MERGED_REF, RS_MERGED_REF, SAMPLING_MERGED_REF, PATTERN_MERGED_REF = my_params[
            f"merged_ref"]
        KMER_SIZE_MERGED_SAMPLE, RS_MERGED_SAMPLE, SAMPLING_MERGED_SAMPLE, PATTERN_MERGED_SAMPLE = my_params[
            f"merged_sample"]
    # if any error happens
    except:
        my_output_msg(format_exc())
        raise ValueError(
            "Incorrect or missing parameters file ; check path and/or contents of json reference.")

    # init iterables and memory spaces
    topmost: dict = {}
    test_results: dict = {}
    Path(f"output/{JOB}/").mkdir(parents=True, exist_ok=True)

    targeted_taxas: list[str] = TAXAS_LEVELS[:TAXAS_LEVELS.index(
        targeted_level)+1]
    output: dict = {level: None for level in targeted_taxas}

    for taxa in targeted_taxas:
        KMER_SIZE_REF, RS_REF, SAMPLING_REF, PATTERN_REF = my_params[f"{taxa}_ref"]
        KMER_SIZE_SAMPLE, RS_SAMPLE, SAMPLING_SAMPLE, PATTERN_SAMPLE = my_params[
            f"{taxa}_sample"]

        list_parent_level = output[
            f"Possible for {TAXAS_LEVELS[TAXAS_LEVELS.index(taxa)-1]}"] if taxa != 'domain' else [False]

        for parent_level in list_parent_level:
            if isinstance(parent_level, bool):
                parent_level = None

############################################ DATABASE STUFF ###############################################

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
                    read_size=RS_REF,
                    classif_level=taxa,
                    sp_determied=parent_level,
                    pattern=PATTERN_REF
                )

            map_sp = load_mapping(OUTPUT_PATH, DATABASE,
                                  taxa, parent_level)

            output[f"{parent_level} diversity"] = list(map_sp.keys())

############################################ MODEL STUFF ###############################################

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
                read_size=RS_SAMPLE,
                classif_level=taxa,
                sp_determied=parent_level,
                pattern=PATTERN_SAMPLE
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

############################################ TEST STUFF ###############################################

            # base tests for heatmap and evaluators
            test_results[f"{taxa}_{parent_level}"] = (test_model(
                OUTPUT_PATH, JOB, DATABASE, taxa, reads_threshold, parent_level, func_reads))

            if single_way:

                output_temp = test_unk_sample(
                    OUTPUT_PATH, JOB, DATABASE, taxa, parent_level, threshold, reads_threshold, test_state, SAMPLING_SAMPLE, func_reads)
                topmost[f"{taxa}_{parent_level}"] = output_temp[f"Reads summation {taxa}"]

            else:
                # add second call for two ways with reverse comp
                raise ValueError("Double pass not implemented yet.")

            if f"Possible for {taxa}" in output:
                output[f"Possible for {taxa}"] = [
                    *output[f"Possible for {taxa}"], *output_temp[f"Possible for {taxa}"]]

            output = {**output_temp, **output}

############################################ MERGED STUFF ###############################################

    if not check_if_merged_database_exists(DATABASE, OUTPUT_PATH):

        make_datasets(
            input_style=False,
            job_name=JOB,
            input_dir=INPUT_PATH,
            path=OUTPUT_PATH,
            datas=['train', 'test'],
            db_name=DATABASE,
            sampling=SAMPLING_MERGED_REF,
            kmer_size=KMER_SIZE_MERGED_REF,
            read_size=RS_MERGED_REF,
            classif_level='merged',
            sp_determied='merged',
            pattern=PATTERN_MERGED_REF
        )

    map_merged_sp = load_mapping(OUTPUT_PATH, DATABASE,
                                 'merged', 'merged')

    if force_rebuild or not check_if_merged_model_exists(DATABASE, OUTPUT_PATH):
        # needs to create merged database
        make_model(JOB, OUTPUT_PATH, 'merged', DATABASE,
                   'merged', init_parameters(len(map_merged_sp), tree_depth), number_rounds=nr)

    # then test in merged mode the sample (sample against all families without splitting steps)
    make_datasets(
        input_style=input_file,
        job_name=JOB,
        input_dir=INPUT_PATH,
        path=OUTPUT_PATH,
        datas=['unk'],
        db_name=DATABASE,
        sampling=SAMPLING_MERGED_SAMPLE,
        kmer_size=KMER_SIZE_MERGED_SAMPLE,
        read_size=RS_MERGED_SAMPLE,
        classif_level='merged',
        sp_determied='merged',
        pattern=PATTERN_MERGED_SAMPLE
    )
    output_merged_sample = test_unk_sample(
        OUTPUT_PATH, JOB, DATABASE, 'merged', 'merged', threshold, reads_threshold, test_state, SAMPLING_MERGED_SAMPLE, func_reads)

    plot_pie_merge(output_merged_sample, JOB)
############################################ END ###############################################

    for i, k in enumerate(targeted_taxas):
        responses: list[str] = [f"domain_None", f"phylum_{output['domain']}",
                                f"group_{output['phylum']}", f"order_{output['group']}", f"family_{output['order']}"]
        output[k] = topmost[responses[i]]

    # extraction of most probable path
    path_taxa = [f"{output[k]} ({k[0]})" for k in targeted_taxas]

    output["parcimonious_path"] = tree_render(output, JOB, path_taxa)

    save_output({'Date': f"{datetime.today().strftime('%Y.%m.%d - %H:%M:%S')}", **
                 vars(args), **output}, JOB)
    #gen_html_report(my_params, JOB, [], output, targeted_taxas,test_results, threshold, test_state, output_merged_sample, round(reads_threshold, 2))
    make_doc(JOB, my_params, TAXAS_LEVELS, output, test_results,
             test_state, threshold, reads_threshold)
