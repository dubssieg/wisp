"Functions to manipulate data such as JSON, LIBSVM..."

from xgboost import DMatrix
from os import listdir
from pathlib import Path
from os.path import exists
from json import load
from python_tools import my_output_msg


def load_mapping(path: str, db_name: str, classif_level: str, sp_determined: str | None) -> dict:
    """Check the classification level we're working at, and load corresponding mapping file

    Args:
        path (str): output path which is in the json parameters file
        db_name (str): name of database user choosed
        classif_level (str): classif level we're currently trying to work on
        sp_determined (str | None): upper taxa if already determined

    Returns:
        dict: mapping from clade names to int
    """
    if sp_determined == None:
        my_path = f"{path}{db_name}/{classif_level}/saved_mapping.json"
    else:
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_saved_mapping.json"
    with open(my_path, "r") as file_manager:
        return load(file_manager)


def species_list(input_dir: str, filter: str | None = None) -> list:
    """Return the list of all species in a given dir, based upon their filenames

    Args:
        input_dir (str): directory were we're checking for genomes to load
        filter (str | None, optional): perviously determined taxa to filter out results. Defaults to None.

    Returns:
        list: all species matching our filter
    """
    if filter == None:
        list_match = [file.split('.')[0] for file in listdir(input_dir)]
        return list_match
    else:
        list_match = [file.split('.')[0] for file in listdir(
            input_dir) if file.find(filter) != -1]
        return list_match


def species_map(input_dir: str, int_level: int, filter: str | None = None) -> dict:
    """Return a dict, association of a species name and a code for species

    Args:
        input_dir (str): _description_
        filter (str | None, optional): _description_. Defaults to None.

    Returns:
        dict: _description_
    """
    species_listed = list(set([file.split('_')[int_level]
                          for file in species_list(input_dir, filter)]))
    return {species_listed[i]: i for i in range(len(species_listed))}


def load_xgboost_data(path: str, classif_level: str, suffix: str, db_name: str, sp_determined: str | None, sample_name: str):
    """
    Return the target dataframe

    * path (str) : location where .txt.train and .txt.test are stored
    * classif_level (str) : level of cassification we're working at
    * suffix (str) : 'train', 'test', 'unk'
    * db_name (str) : database we need to search in
    * sp_determined (str | None): upper level we've already determined
    """
    if suffix == 'unk':
        my_path = f"{path}data.txt.{suffix}"
    elif sp_determined == None:
        my_path = f"{path}{db_name}/{classif_level}/data.txt.{suffix}"
    else:
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_data.txt.{suffix}"
    return DMatrix(my_path)


def write_xgboost_data(data: list[str], path: str, classif_level: str, suffix: str, db_name: str, sample_name: str, sp_determined: str | None) -> None:
    """
    Writes out the xgboost data

    * data (list[str]) : a list of strings in LIBSVM format
    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * suffix (str) : 'train', 'test', 'unk'
    * db_name (str) : database we need to search in
    * sp_determined (str | None): upper level we've already determined
    """
    if suffix == 'unk':
        Path(f"{path}").mkdir(parents=True, exist_ok=True)
        my_path = f"{path}data.txt.{suffix}"
    elif sp_determined == None:
        Path(f"{path}{db_name}/{classif_level}/").mkdir(parents=True, exist_ok=True)
        my_path = f"{path}{db_name}/{classif_level}/data.txt.{suffix}"
    else:
        Path(f"{path}{db_name}/{classif_level}/").mkdir(parents=True, exist_ok=True)
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_data.txt.{suffix}"
    with open(my_path, "w") as writer:
        writer.write('\n'.join(data))


def check_if_merged_database_exists(db_name: str, path: str) -> bool:
    """Checks if train and test files exists for current situation

    Args:
        db_name (str): name of db, user imput
        path (str): path to the database, user input from parameters json file
        taxa_level (str): level we're working at
        sp_determined (str | None): upper level we've already determined

    Returns:
        bool: if all files are present in targeted folder
    """
    files_to_check = [f"{path}{db_name}/merged/merged_data.txt.train",
                      f"{path}{db_name}/merged/merged_data.txt.test", f"{path}{db_name}/merged/merged_saved_mapping.json"]
    all_correct = True
    for file in files_to_check:
        if not exists(file):
            my_output_msg(f"File not existing : {file}")
            all_correct = False
        else:
            my_output_msg(f"File found : {file}")
    return all_correct


def check_if_database_exists(db_name: str, path: str, taxa_level: str, sp_determined: str | None) -> bool:
    """Checks if train and test files exists for current situation

    Args:
        db_name (str): name of db, user imput
        path (str): path to the database, user input from parameters json file
        taxa_level (str): level we're working at
        sp_determined (str | None): upper level we've already determined

    Returns:
        bool: if all files are present in targeted folder
    """
    files_to_check = [f"{path}{db_name}/{taxa_level}/data.txt.train", f"{path}{db_name}/{taxa_level}/data.txt.test", f"{path}{db_name}/{taxa_level}/saved_mapping.json"] if sp_determined == None else [
        f"{path}{db_name}/{taxa_level}/{sp_determined}_data.txt.train", f"{path}{db_name}/{taxa_level}/{sp_determined}_data.txt.test", f"{path}{db_name}/{taxa_level}/{sp_determined}_saved_mapping.json"]
    all_correct = True
    for file in files_to_check:
        if not exists(file):
            my_output_msg(f"File not existing : {file}")
            all_correct = False
        else:
            my_output_msg(f"File found : {file}")
    return all_correct


def check_if_model_exists(db_name: str, path: str, taxa_level: str, sp_determined: str | None) -> bool:
    """Checks if a XGBoost model as already been build for current situation

    Args:
        db_name (str): name of db, user imput
        path (str): path to the database, user input from parameters json file
        taxa_level (str): level we're working at
        sp_determined (str | None): upper level we've already determined

    Returns:
        bool: if all files are present in targeted folder
    """
    files_to_check = [f"{path}{db_name}/{taxa_level}/saved_params.json", f"{path}{db_name}/{taxa_level}/saved_model.json"] if sp_determined == None else [
        f"{path}{db_name}/{taxa_level}/{sp_determined}_saved_params.json", f"{path}{db_name}/{taxa_level}/{sp_determined}_saved_model.json"]
    all_correct = True
    for file in files_to_check:
        if not exists(file):
            my_output_msg(f"File not existing : {file}")
            all_correct = False
        else:
            my_output_msg(f"File found : {file}")
    return all_correct


def check_if_merged_model_exists(db_name: str, path: str) -> bool:
    """Checks if a XGBoost merged model (whole family collection) as already been build for current situation

    Args:
        db_name (str): name of db, user imput
        path (str): path to the database, user input from parameters json file
        taxa_level (str): level we're working at
        sp_determined (str | None): upper level we've already determined

    Returns:
        bool: if all files are present in targeted folder
    """
    files_to_check = [f"{path}{db_name}/merged/merged_saved_params.json",
                      f"{path}{db_name}/merged/merged_saved_model.json"]
    all_correct = True
    for file in files_to_check:
        if not exists(file):
            my_output_msg(f"File not existing : {file}")
            all_correct = False
        else:
            my_output_msg(f"File found : {file}")
    return all_correct
