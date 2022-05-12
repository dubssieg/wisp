import argparse
from sample_class import make_datasets
from build_softprob import load_model, load_params, prediction, make_model, init_parameters, make_testing
import warnings
import xgboost as xgb
import json
from python_tools import my_output_msg, my_logs_global_config
from collections import Counter
from datetime import datetime
from wisp_view import reads_species_plotter, gen_html_report, compare_test, tree_render, plot_boosting
from wisp_lib import load_xgboost_data, check_if_database_exists, check_if_model_exists, load_mapping, load_json
from pathlib import Path


def estimations(preds, sample_name: str, inverted_map: dict[str, str], clade: str, determined: str, threshold: float = 0.1) -> dict:
    """
    Does calculation upon number of predicated reads
    Returns a dict and call for plotting results as a barplot

    * preds (numpy.array) : contains all teh preds made by the model
    * sample_name (str) : name of job, used to name plots
    * inverted_map (dict) : mapping numbers of classes to taxas
    * clade (str) : level of taxonmy we're looking at
    * threshold (float) : limits to percentage subclasses we will analyse
    """
    sum_preds, sum = Counter(
        [int(round(predict, 0)) for predict in preds]), len(preds)

    outputs_predictions: dict = {
        f"{clade} {inverted_map[str(k)]}": f"{v} reads out of {sum} ({round(v/sum*100,2)}%)" for k, v in sum_preds.items() if str(k) in inverted_map.keys()}

    outputs_labels: dict = {
        f"{inverted_map[str(k)]} ({clade[0]})": f"{round(v/sum*100,2)}%" for k, v in sum_preds.items() if str(k) in inverted_map.keys()}

    outputs_classif = {
        f"Reads summation {clade}": inverted_map[str(max(sum_preds, key=sum_preds.get))]
    }

    map_clade = ['-', 'd', 'p', 'g', 'o', 'f']

    outputs_classif[f"Possible for {taxa}"] = [inverted_map[str(
        k)] for k, v in sum_preds.items() if float(v) > float(threshold * sum)]
    outputs_classif[f"Tree {determined} ({map_clade[map_clade.index(clade[0])-1]})"] = [
        f"{inverted_map[str(k)]} ({clade[0]})" for k, v in sum_preds.items() if float(v) > float(threshold * sum)]

    outputs = {**outputs_classif, **outputs_predictions, **outputs_labels}

    Path(f"output/{JOB}").mkdir(parents=True, exist_ok=True)
    reads_species_plotter(preds, sample_name, inverted_map,
                          clade, determined, threshold)

    return outputs


def save_output(dico: dict, job_name: str) -> None:
    """
    Saves a dict as a .json file

    * dico (dict) : dictionnary to store
    * job_name (str) : name of current job, creates output path
    """
    with open(f"output/{job_name}/{job_name}_results.json", "a") as fm:
        json.dump(dico, fm)
        fm.write('\n')


def test_model(out_path, job_name, database_name, classif_level, reads_threshold, sp_determined: str | None):
    """
    Does the testing of our model with the data in /test.
    Will estimate some meaningful estimators and will plot heatmap

    * out_path (str) : path where to store graphs
    * job_name (str) : name of current job, will give name to files
    * database_name (str) : path to train and test indexes
    * classif_level (str) : current taxa level
    * sp_determined (str) : estimation which was done at taxa t-1
    """
    map_sp = load_mapping(out_path, database_name,
                          classif_level, sp_determined)
    inverted_map = {str(v): k for k, v in map_sp.items()}

    bst = xgb.Booster()

    my_output_msg("Data loading...")
    dtest = load_xgboost_data(
        out_path, classif_level, 'test', database_name, sp_determined)

    my_output_msg("Model loading...")
    load_model(bst, out_path, classif_level, database_name, sp_determined)
    load_params(bst, out_path, classif_level, database_name, sp_determined)

    my_output_msg("Preds calculation...")
    preds = prediction(dtest, bst, job_name, classif_level,
                       sp_determined, reads_threshold, False, inverted_map)

    if sp_determined == None:
        with open(f"{out_path}{database_name}/{classif_level}/data.txt.test", "r") as reader:
            real = [int(l.split(' ')[0]) for l in reader]
    else:
        with open(f"{out_path}{database_name}/{classif_level}/{sp_determined}_data.txt.test", "r") as reader:
            real = [int(l.split(' ')[0]) for l in reader]

    # filtering
    real = [p for i, p in enumerate(real) if not isinstance(preds[i], bool)]
    preds = [p for p in preds if not isinstance(p, bool)]

    return compare_test(real, preds, inverted_map, job_name, classif_level, sp_determined)


def test_unk_sample(out_path, job_name, database_name, classif_level, sp_determined, threshold, reads_threshold, test_status):
    map_sp = load_mapping(out_path, database_name,
                          classif_level, sp_determined)
    inverted_map = {str(v): k for k, v in map_sp.items()}

    bst = xgb.Booster()

    my_output_msg("Data loading...")
    dunk = load_xgboost_data(
        out_path, classif_level, 'unk', database_name, sp_determined)

    my_output_msg("Model loading...")
    load_model(bst, out_path, classif_level, database_name, sp_determined)
    load_params(bst, out_path, classif_level, database_name, sp_determined)

    my_output_msg("Preds calculation...")
    preds = prediction(dunk, bst, job_name, classif_level,
                       sp_determined, reads_threshold, test_status, inverted_map)
    preds = [p for p in preds if not isinstance(p, bool)]

    return estimations(preds, job_name, inverted_map, classif_level, sp_determined, threshold)


if __name__ == "__main__":
    "Executes main procedure"
    warnings.filterwarnings('ignore')  # to ignore xgboost warnnings
    my_logs_global_config("LOG_wisp")

    parser = argparse.ArgumentParser()

    # declaring args
    # parser.add_argument("genome_path", help="genome/read to evaluate", type=str)
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
    # if any error happens
    except:
        raise ValueError("Incorrect or missing file")

    FUNC: str | None = None
    RATIO: float = 1.5

    # init iterables and memory spaces
    topmost: dict = {}
    test_results: dict = {}
    taxa_map: list = ['domain', 'phylum', 'group', 'order', 'family']
    output: dict = {'domain': None, 'phylum': None,
                    'group': None, 'order': None, 'family': None}

    # if only 4 levels, depreciated
    #taxa_map: list = ['domain', 'phylum', 'order', 'family']
    #output: dict = {'domain': None, 'phylum': None,'order': None, 'family': None}
    for taxa in taxa_map:
        KMER_SIZE_REF, RS_REF, SAMPLING_REF = my_params[f"{taxa}_ref"]
        KMER_SIZE_SAMPLE, RS_SAMPLE, SAMPLING_SAMPLE = my_params[f"{taxa}_sample"]

        list_parent_level = output[
            f"Possible for {taxa_map[taxa_map.index(taxa)-1]}"] if taxa != 'domain' else [False]

        for parent_level in list_parent_level:
            if isinstance(parent_level, bool):
                parent_level = None

            # parent_level: str | None = output[taxa_map[taxa_map.index(taxa)-1]] if taxa != 'domain' else None

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
                           parent_level, init_parameters(len(map_sp)), number_rounds=nr)

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

    tree_render(output, JOB, path_taxa)

    save_output({'Date': f"{datetime.today().strftime('%Y.%m.%d - %H:%M:%S')}", **
                 vars(args), **output}, JOB)
    gen_html_report(my_params, JOB, [], output, taxa_map,
                    test_results, threshold, test_state, round(reads_threshold, 2))
