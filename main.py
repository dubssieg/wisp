"Executes main procedure"
from argparse import ArgumentParser
from traceback import format_exc
from sample_class import make_datasets, make_unk_datasets
from build_softprob import make_model, init_parameters, make_testing
from warnings import filterwarnings
from python_tools import my_output_msg, my_fasta_parser
from datetime import datetime
from wisp_view import plot_boosting, make_doc, global_sample_report
from wisp_lib import counter_ultrafast, kmer_indexing_canonical, reverse_comp, splitting_generator, check_if_database_exists, check_if_model_exists, check_if_merged_database_exists, load_mapping, load_json, check_if_merged_model_exists
from predictors import test_unk_sample, save_output, test_model
from pathlib import Path
from _version import get_versions
__version__ = get_versions()['version']
del get_versions

if __name__ == "__main__":
    ############################################ ARGPARSE ###############################################

    filterwarnings('ignore')  # to ignore xgboost warnnings

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

    # we try to load params file and gather data from it
    try:
        my_params: dict = load_json(args.params)
        # storing args

        DATABASE: str = args.database_name
        TAXAS_LEVELS: list[str] = my_params['levels_list']
        TRAIN_PATH: str = my_params['input_train']
        UNK_PATH: str = my_params['input_unk']
        DATABASE_PATH: str = my_params['database_output']
        REPORTS_PATH: str = my_params['reports_output']
        WINDOW: int = my_params['window_size']
        SAMPLING_OBJECTIVE: int = my_params['sampling_objective']

        nr = int(my_params['nb_boosts'])
        threshold = float(my_params['threshold'])
        reads_threshold = float(my_params['reads_th'])
        test_state = str(my_params['test_mode'])
        force_rebuild = bool(my_params['force_model_rebuild'])
        tree_depth = int(my_params['tree_depth'])
        func_reads = str(my_params['selection_mode'])
        single_way = bool(my_params['single_way'])
        targeted_level = str(my_params['targeted_level'])
        KMER_SIZE_MERGED_REF, SAMPLING_MERGED_REF, PATTERN_MERGED_REF = my_params[
            f"merged_ref"]
        KMER_SIZE_MERGED_SAMPLE, PATTERN_MERGED_SAMPLE = my_params[
            f"merged_sample"]
    # if any error happens
    except Exception as exc:
        my_output_msg(format_exc())
        raise ValueError(
            "Incorrect or missing parameters file ; check path and/or contents of json reference.") from exc

    input_file: str = f"{args.file}"
    fasta_reads: dict = my_fasta_parser(
        f"{UNK_PATH}{input_file}")

    all_reads_report: dict = {}  # will contain all reports, read by read
    global_reports_path: str = f"{REPORTS_PATH}{input_file}/"
    Path(global_reports_path).mkdir(parents=True, exist_ok=True)
    for i, (id, read) in enumerate(fasta_reads.items()):
        if len(read) < WINDOW:
            raise ValueError(
                "Read size is too short given the current window size ; skipping...")
        else:
            # read is correct & clean, and we have good params
            if single_way:
                my_output_msg(
                    f"Processing read {i+1} out of {len(fasta_reads)}...")
                jobs: list[str] = [f"{args.job_name}_read_{i}"]
                path_for_read: list[str] = [
                    f"{REPORTS_PATH}{input_file}/{job}/" for job in jobs]
                all_reads = [read for read in splitting_generator(
                    read, WINDOW, SAMPLING_OBJECTIVE)]
                count_func = counter_ultrafast
            else:
                my_output_msg(
                    f"Processing read {i+1} out of {len(fasta_reads)}...")
                jobs: list[str] = [
                    f"{args.job_name}_read_{i}_strand_A", f"{args.job_name}_read_{i}_strand_B"]
                path_for_read: list[str] = [
                    f"{REPORTS_PATH}{input_file}/{job}/" for job in jobs]
                rds = [read, reverse_comp(read)]
                all_reads = [splitting_generator(
                    read, WINDOW, SAMPLING_OBJECTIVE) for read in rds]
                count_func = counter_ultrafast
            for path in path_for_read:
                Path(path).mkdir(parents=True, exist_ok=True)

            # init iterables and memory spaces
            topmost: dict = {}
            test_results: dict = {}

            targeted_taxas: list[str] = TAXAS_LEVELS[:TAXAS_LEVELS.index(
                targeted_level)+1]
            output: dict = {level: None for level in targeted_taxas}

            for i, JOB in enumerate(jobs):

                for taxa in targeted_taxas:
                    KMER_SIZE_REF, SAMPLING_REF, PATTERN_REF = my_params[f"{taxa}_ref"]
                    KMER_SIZE_SAMPLE, PATTERN_SAMPLE = my_params[f"{taxa}_sample"]

                    list_parent_level = output[
                        f"Possible for {TAXAS_LEVELS[TAXAS_LEVELS.index(taxa)-1]}"] if taxa != 'domain' else [False]

                    for parent_level in list_parent_level:
                        if isinstance(parent_level, bool):
                            parent_level = None

            ############################################ DATABASE STUFF ###############################################

                        if not check_if_database_exists(DATABASE, DATABASE_PATH, taxa, parent_level):

                            make_datasets(
                                input_style=False,
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

                        map_sp = load_mapping(DATABASE_PATH, DATABASE,
                                              taxa, parent_level)

                        output[f"{parent_level} diversity"] = list(
                            map_sp.keys())

            ############################################ MODEL STUFF ###############################################

                        if force_rebuild or not check_if_model_exists(DATABASE, DATABASE_PATH, taxa, parent_level):

                            make_model(JOB, DATABASE_PATH, taxa, DATABASE,
                                       parent_level, init_parameters(len(map_sp), tree_depth), number_rounds=nr)

                        number_of_reads = make_unk_datasets(
                            func=count_func,
                            all_reads=all_reads,
                            job_name=JOB,
                            path=DATABASE_PATH,
                            db_name=DATABASE,
                            kmer_size=KMER_SIZE_SAMPLE,
                            classif_level=taxa,
                            sp_determied=parent_level,
                            pattern=PATTERN_SAMPLE
                        )

                        # full test set, takes time, but gives info on structure
                        if test_state == 'verbose':
                            successive_boost_results = make_testing(
                                path_to_save=path_for_read[i],
                                size_kmer=KMER_SIZE_REF,
                                job_name=JOB,
                                sp_determined=parent_level,
                                path=DATABASE_PATH,
                                db_name=DATABASE,
                                classif_level=taxa,
                                class_count=len(map_sp),
                                model_parameters=init_parameters(len(map_sp)),
                                number_rounds=nr
                            )

                            plot_boosting(successive_boost_results,
                                          JOB, taxa, parent_level, nr)

            ############################################ TEST STUFF ###############################################
                        if test_state != 'no_test':
                            # base tests for heatmap and evaluators
                            test_results[f"{taxa}_{parent_level}"] = (test_model(path_for_read[i],
                                                                                 DATABASE_PATH, JOB, DATABASE, taxa, reads_threshold, parent_level, func_reads))

                        output_temp = test_unk_sample(path_for_read[i],
                                                      DATABASE_PATH, JOB, DATABASE, taxa, parent_level, threshold, reads_threshold, test_state, number_of_reads, func_reads, test_state)
                        topmost[f"{taxa}_{parent_level}"] = output_temp[f"Reads summation {taxa}"]

                        if f"Possible for {taxa}" in output:
                            output[f"Possible for {taxa}"] = [
                                *output[f"Possible for {taxa}"], *output_temp[f"Possible for {taxa}"]]

                        output = {**output_temp, **output}

            ############################################ MERGED STUFF ###############################################
                """ Currently bugged
                if not check_if_merged_database_exists(DATABASE, DATABASE_PATH):

                    make_datasets(
                        input_style=False,
                        job_name=JOB,
                        input_dir=TRAIN_PATH,
                        path=DATABASE_PATH,
                        datas=['train', 'test'],
                        db_name=DATABASE,
                        sampling=SAMPLING_MERGED_REF,
                        kmer_size=KMER_SIZE_MERGED_REF,
                        read_size=WINDOW,
                        classif_level='merged',
                        sp_determied='merged',
                        pattern=PATTERN_MERGED_REF
                    )

                map_merged_sp = load_mapping(DATABASE_PATH, DATABASE,
                                             'merged', 'merged')

                if force_rebuild or not check_if_merged_model_exists(DATABASE, DATABASE_PATH):
                    # needs to create merged database
                    make_model(JOB, DATABASE_PATH, 'merged', DATABASE,
                               'merged', init_parameters(len(map_merged_sp), tree_depth), number_rounds=nr)

                # then test in merged mode the sample (sample against all families without splitting steps)
                make_unk_datasets(
                    func=count_func,
                    all_reads=all_reads[i],
                    job_name=JOB,
                    path=DATABASE_PATH,
                    db_name=DATABASE,
                    kmer_size=KMER_SIZE_MERGED_SAMPLE,
                    classif_level='merged',
                    sp_determied='merged',
                    pattern=PATTERN_MERGED_SAMPLE
                )

                if test_state != 'no_test':
                    output_merged_sample = test_unk_sample(path_for_read[i],
                                                           DATABASE_PATH, JOB, DATABASE, 'merged', 'merged', threshold, reads_threshold, test_state, len(all_reads[i]), func_reads, test_state)
                    plot_pie_merge(path_for_read[i], output_merged_sample, JOB)
                """
            ############################################ END ###############################################

                for ib, k in enumerate(targeted_taxas):
                    responses: list[str] = [f"domain_None", f"phylum_{output['domain']}",
                                            f"group_{output['phylum']}", f"order_{output['group']}", f"family_{output['order']}"]
                    output[k] = topmost[responses[ib]]

                # extraction of most probable path
                path_taxa = [f"{output[k]} ({k[0]})" for k in targeted_taxas]

                #output["parcimonious_path"] = tree_render(path_for_read, output, JOB, path_taxa)

                all_reads_report[id] = [output[k] for k in targeted_taxas]

                save_output({'Date': f"{datetime.today().strftime('%Y.%m.%d - %H:%M:%S')}", **
                            vars(args), **output}, JOB, path_for_read[i])
                if test_state != 'no_test':
                    mts = True if test_state == 'verbose' else False
                    make_doc(path_for_read[i], JOB, my_params, TAXAS_LEVELS, output, test_results,
                             mts, threshold, reads_threshold, len(all_reads[i]), id, __version__, len(read))

    # saves global output (merged reads attribution)
    save_output({'Date': f"{datetime.today().strftime('%Y.%m.%d - %H:%M:%S')}",
                'Job': input_file, **all_reads_report}, input_file, global_reports_path)
    global_sample_report(global_reports_path, input_file, all_reads_report)

    with open(f"{global_reports_path}{input_file}.txt", 'w') as writer:
        writer.write(
            '\n'.join([f"{k}:{v}" for k, v in all_reads_report.items()]))
