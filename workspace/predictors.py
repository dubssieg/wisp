from build_softprob import load_model, load_params, prediction
import xgboost as xgb
from json import dump
from wisp_tools import my_output_msg
from collections import Counter
from wisp_view import reads_species_plotter, compare_test
from wisp_lib import load_xgboost_data, load_mapping


def estimations_merged(preds, sample_name: str, inverted_map: dict[str, str], clade: str, determined: str, threshold: float, sampling_number: int) -> dict:
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
        f"{inverted_map[str(k)]}": v/sum for k, v in sum_preds.items() if str(k) in inverted_map.keys()}

    return {**outputs_predictions}


def estimations(path_to_save, preds, sample_name: str, inverted_map: dict[str, str], clade: str, determined: str, threshold: float, sampling_number: int, test_state: str) -> dict:
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

    outputs_classif[f"Possible for {clade}"] = [inverted_map[str(
        k)] for k, v in sum_preds.items() if float(v) > float(threshold * sum)]
    outputs_classif[f"Tree {determined} ({map_clade[map_clade.index(clade[0])-1]})"] = [
        f"{inverted_map[str(k)]} ({clade[0]})" for k, v in sum_preds.items() if float(v) > float(threshold * sum)]
    if test_state == 'min_set' or test_state == 'verbose':
        reads_species_plotter(path_to_save, preds, sample_name, inverted_map,
                              clade, determined, threshold, sampling_number)

    return {**outputs_classif, **outputs_predictions, **outputs_labels}


def save_output(dico: dict, job_name: str, path_for_read: str) -> None:
    """
    Saves a dict as a .json file

    * dico (dict) : dictionnary to store
    * job_name (str) : name of current job, creates output path
    """
    with open(f"{path_for_read}{job_name}_results.json", "a") as fm:
        dump(dico, fm)
        fm.write('\n')


def test_model(path_to_save: str, out_path, job_name, database_name, classif_level, reads_threshold, sp_determined: str | None, func, excluded):
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
        out_path, classif_level, 'test', database_name, sp_determined, job_name, exclude=excluded)

    my_output_msg("Model loading...")
    load_model(bst, out_path, classif_level,
               database_name, sp_determined, job_name)
    load_params(bst, out_path, classif_level, database_name, sp_determined)

    my_output_msg("Preds calculation...")
    preds = prediction(dtest, bst, job_name, classif_level,
                       sp_determined, reads_threshold, False, inverted_map, func, True)

    if excluded != []:
        print(excluded)
        test_path = f"{out_path}{database_name}/temp/{job_name}/data.txt.test"
    elif sp_determined is None:
        test_path = f"{out_path}{database_name}/{classif_level}/data.txt.test"
    else:
        test_path = f"{out_path}{database_name}/{classif_level}/{sp_determined}_data.txt.test"
    with open(test_path, "r") as reader:
        real = [int(l.split(' ')[0])
                for l in reader if l.split(' ')[0] != '\n']

    # filtering
    real = [p for i, p in enumerate(real) if not isinstance(preds[i], bool)]
    preds = [p for p in preds if not isinstance(p, bool)]

    return compare_test(path_to_save, real, preds, inverted_map, job_name, classif_level, sp_determined)


def test_unk_sample(path_to_save, out_path, job_name, database_name, classif_level, sp_determined, threshold, reads_threshold, test_status, sampling_number, func, test_state, bool_temporary):
    map_sp = load_mapping(out_path, database_name,
                          classif_level, sp_determined)
    inverted_map = {str(v): k for k, v in map_sp.items()}

    bst = xgb.Booster()

    my_output_msg("Model loading...")
    load_model(bst, out_path, classif_level, database_name,
               sp_determined, job_name, bool_temporary)
    load_params(bst, out_path, classif_level, database_name, sp_determined)

    # loop over all data to exit model re-loading at each step?

    my_output_msg("Data loading...")
    dunk = load_xgboost_data(
        out_path, classif_level, 'unk', database_name, sp_determined, job_name)

    my_output_msg("Preds calculation...")
    preds = prediction(dunk, bst, job_name, classif_level,
                       sp_determined, reads_threshold, test_status, inverted_map, func, True)
    preds = [p for p in preds if not isinstance(p, bool)]

    if classif_level != 'merged':
        return estimations(path_to_save, preds, job_name, inverted_map, classif_level, sp_determined, threshold, sampling_number, test_state)
    else:
        return estimations_merged(preds, job_name, inverted_map, classif_level, sp_determined, threshold, sampling_number)

    # end loop
