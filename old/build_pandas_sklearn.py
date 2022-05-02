import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import xgboost as xgb
from python_tools import my_output_msg, my_function_timer
from wisp_lib import load_xgboost_data
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split

# consts for test
DEPTH: int = 6
ETA: float = 0.3
TREES: int = 6


def plot_confusion_matrix(cm, classes, normalized=True, cmap='bone'):
    plt.figure(figsize=(7, 6))
    norm_cm = cm
    if normalized:
        norm_cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        sns.heatmap(norm_cm, annot=cm, fmt='g', xticklabels=classes,
                    yticklabels=classes, cmap=cmap)

#################### INIT FUNCTIONS ####################


def global_scope():
    xgb.config_context(booster='gbtree', min_child_weight=1,
                       tree_method='approx', predictor='cpu_predictor')


def init_parameters(depth: int = DEPTH, eta: float = ETA):
    """
    Store algorithm args into dict in order to control our XGBoost model.
    Return a dict of args 

    * depth (int)
    * eta (int)
    """
    return {'max_depth': depth, 'eta': eta}

#################### SAVE AND LOAD ####################


def save_model(bst, path: str, classif_level: str, db_name: str, sp_determined):
    """
    Store model in order to dodge calculation at each new run
    Executes a memory snapshot and put it in a .json file

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * bst (xgb.Model) : xgboost model
    * db_name (str) : database we need to search in
    """
    if sp_determined == None:
        my_path = f"{path}{db_name}/{classif_level}/saved_model.json"
    else:
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_saved_model.json"
    bst.save_model(my_path)


def load_model(bst, path: str, classif_level: str, db_name: str, sp_determined):
    """
    Load previously stored model in order to dodge calculation
    Load a memory snapshot from a .json file

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * bst (xgb.Model) : xgboost model
    * db_name (str) : database we need to search in
    """
    if sp_determined == None:
        my_path = f"{path}{db_name}/{classif_level}/saved_model.json"
    else:
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_saved_model.json"
    bst.load_model(my_path)


def save_params(bst, path: str, classif_level: str, db_name: str, sp_determined):
    """
    Store params of a model in order to dodge configuration on other model

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * bst (xgb.Model) : xgboost model
    * db_name (str) : database we need to search in
    """
    if sp_determined == None:
        my_path = f"{path}{db_name}/{classif_level}/saved_params.json"
    else:
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_saved_params.json"
    with open(my_path, "w") as writer:
        writer.write(bst.save_config())


def load_params(bst, path: str, classif_level: str, db_name: str, sp_determined):
    """
    Load params of a model in order to dodge configuration on other model

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * bst (xgb.Model) : xgboost model
    * db_name (str) : database we need to search in
    """
    if sp_determined == None:
        my_path = f"{path}{db_name}/{classif_level}/saved_params.json"
    else:
        my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_saved_params.json"
    with open(my_path, "r") as reader:
        bst.load_config(''.join([line for line in reader]))

#################### TRAINING ####################


def modelisation(dtrain, params: dict, num_round: int):
    """
    Does the training of our model (here to test)
    Return trained model

    * dtrain (DMatrix) : data on which model will be established
    * params (dict) : parameters to use for the model
    * num_round : number of successives boostings (default value 2)
    """
    return xgb.train(params, dtrain, num_round)


def prediction(data, model):
    """
    Does the calculation given a model and a dataset of catergories of this dataset
    Return predictions

    * data (DMatrix) : data to compute
    * model (xgb.Booster) : trained XGBoost model
    """
    return model.predict(data)

#################### CORE ####################


@my_function_timer("Building model")
def make_model(path: str, classif_level: str, db_name: str, sp_determined, model_parameters: dict = init_parameters(), number_rounds: int = 10):
    """
    Main procedure : builds the model and saves it

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * db_name (str) : database we need to search in   
    * model_parameters (dict) : a dict of parameters for our xgboost model
    * number_rounds (int) : number of rounds we go through to boost our trees 
    """
    bst = xgb.Booster()
    global_scope()
    dtrain = load_xgboost_data(
        path, classif_level, 'train', db_name, sp_determined)
    bst = modelisation(dtrain, model_parameters, number_rounds)
    save_model(bst, path, classif_level, db_name, sp_determined)
    save_params(bst, path, classif_level, db_name, sp_determined)
