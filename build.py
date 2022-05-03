import xgboost as xgb
from python_tools import my_output_msg, my_function_timer
from wisp_lib import load_xgboost_data, recode_kmer_4
from wisp_view import compare_test, plot_features
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt

# consts for test
DEPTH: int = 6
ETA: float = 0.3
TREES: int = 6

#################### INIT FUNCTIONS ####################


def global_scope():
    """
    Set up global xgboost parameters
    """
    xgb.config_context(booster='gbtree', min_child_weight=1,
                       tree_method='approx', predictor='cpu_predictor')


def init_parameters(class_count: int, depth: int = DEPTH, eta: float = ETA):
    """
    Store algorithm args into dict in order to control our XGBoost model.
    Return a dict of args 

    * depth (int)
    * eta (int)
    """
    return {'max_depth': depth, 'objective': 'multi:softmax', 'num_class': class_count, 'eta': eta}

#################### VALIDATION ####################


def tests_with_kfold(cls, xgbMatrix, X, y, class_count, classif_level, inverted_map, sample, sp_determined):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=0)

    # cls = xgb.XGBClassifier(**init_parameters(class_count))

    cls.fit(X_train, y_train)

    test_preds = cls.predict(X_test)

    xgb_cv = xgb.cv(dtrain=xgbMatrix, params=init_parameters(class_count), nfold=3,
                    num_boost_round=50, early_stopping_rounds=10, metrics="auc", as_pandas=True, seed=123)

    xgb.plot_importance(cls)
    plt.figure(figsize=(16, 12))
    plt.savefig(f"{classif_level}_feature_importance.png")

    return (accuracy_score(y_test, test_preds), xgb_cv, compare_test(y_test, test_preds, inverted_map, sample, classif_level, sp_determined))


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
def make_model(job_name: str, path: str, classif_level: str, db_name: str, sp_determined, model_parameters: dict, number_rounds: int):
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
    #X, y = pandas_load_libsvm(sp_determined, 'train',path, db_name, classif_level)
    dtrain = load_xgboost_data(
        path, classif_level, 'train', db_name, sp_determined)
    #dtrain = xgb.DMatrix(data=X, label=y)
    bst = modelisation(dtrain, model_parameters, number_rounds)
    save_model(bst, path, classif_level, db_name, sp_determined)
    save_params(bst, path, classif_level, db_name, sp_determined)


@my_function_timer("Plotting feature importance")
def make_testing(size_kmer, job_name, sp_determined, path, db_name, classif_level, class_count, model_parameters, number_rounds) -> tuple:
    """_summary_

    Args:
        sp_determined (_type_): _description_
        suffix (_type_): _description_
        path (_type_): _description_
        db_name (_type_): _description_
        classif_level (_type_): _description_
        inverted_map (_type_): _description_
        sample (_type_): _description_
        class_count (_type_): _description_

    Returns:
        tuple: tuple containing :
            [0] -> accuracy score
            [1] -> cross-validation
            [2] -> parameters estimations
    """
    bst = xgb.Booster()
    global_scope()
    mat = load_xgboost_data(
        path, classif_level, 'train', db_name, sp_determined)
    xgb_cv = xgb.cv(dtrain=mat, params=init_parameters(class_count), nfold=3,
                    num_boost_round=20, early_stopping_rounds=10, metrics="auc", as_pandas=True, seed=123)
    bst = modelisation(mat, model_parameters, number_rounds)
    mapped = {
        f"{recode_kmer_4(str(k[1:]),size_kmer)}": v for k, v in bst.get_score(importance_type='gain').items()}
    plot_features(mapped, job_name, classif_level, sp_determined)

    """
    tree = xgb.to_graphviz(bst)
    tree.graph_attr = {'dpi': '400'}
    tree.render(
        f"output/{job_name}/{classif_level}_{sp_determined}_trees_overview", format='png')
    """
    return xgb_cv
