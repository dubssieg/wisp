from logging import warning
import xgboost as xgb
from python_tools import my_output_msg, my_function_timer
from wisp_lib import load_xgboost_data, recode_kmer_4
from wisp_view import compare_test, plot_features, plot_all_reads
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from numpy import argmax, amax, mean
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd

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
    return {'max_depth': depth, 'objective': 'multi:softprob', 'num_class': class_count, 'eta': eta, 'eval_metric': 'mlogloss'}

#################### VALIDATION ####################


def tests_with_kfold(path_to_save, cls, xgbMatrix, X, y, class_count, classif_level, inverted_map, sample, sp_determined):
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
    # here for reads division according to species
    return (accuracy_score(y_test, test_preds), xgb_cv, compare_test(path_to_save, y_test, test_preds, inverted_map, sample, classif_level, sp_determined))


#################### SAVE AND LOAD ####################

def save_model(bst, path: str, classif_level: str, db_name: str, sp_determined: str | None, job_name: str, bool_temporary: bool):
    """
    Store model in order to dodge calculation at each new run
    Executes a memory snapshot and put it in a .json file

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * bst (xgb.Model) : xgboost model
    * db_name (str) : database we need to search in
    """
    if bool_temporary:
        my_path = f"{path}{db_name}/temp/{job_name}/saved_model.json"
    else:
        if sp_determined is None:
            my_path = f"{path}{db_name}/{classif_level}/saved_model.json"
        else:
            my_path = f"{path}{db_name}/{classif_level}/{sp_determined}_saved_model.json"
    my_output_msg(f"Model saved in {my_path}", warning)
    bst.save_model(my_path)


def load_model(bst, path: str, classif_level: str, db_name: str, sp_determined: str | None, job_name: str, bool_temporary: bool = False):
    """
    Load previously stored model in order to dodge calculation
    Load a memory snapshot from a .json file

    * path (str) : location where .txt.train and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * bst (xgb.Model) : xgboost model
    * db_name (str) : database we need to search in
    """
    if bool_temporary:
        my_path = f"{path}{db_name}/temp/{job_name}/saved_model.json"
    else:
        if sp_determined is None:
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


def prediction(data, model, sample_name, clade, determined, reads_threshold, test_state, inverted_map, func, with_softmax_norm):
    """
    Does the calculation given a model and a dataset of catergories of this dataset
    Return predictions

    * data (DMatrix) : data to compute
    * model (xgb.Booster) : trained XGBoost model
    """
    preds = model.predict(data)
    res = softmax_from_prediction(preds, reads_threshold, func)
    if test_state == 'verbose':
        softpred_from_prediction(
            preds, sample_name, clade, determined, inverted_map)
        plot_all_reads(preds, sample_name, inverted_map, clade, determined)
    return res if with_softmax_norm else preds


def softmax_from_prediction(preds, reads_selection_threshold, func):
    """Given a list of predictions of x arrays of size num_classes, computes a list of x predictions, one for each sample
    Recursive function that will compute until a selection of reads can be made if threshold is too high

    Args:
        preds (array): a prediction from a objective = softprob model
        reads_selection_threshold (float) : a threshold for meaningful reads selection

    Returns:
        list : a list of x softmax predictions selected from reads above threshold
    """
    if reads_selection_threshold <= 0:
        return [argmax(a) for a in preds]
    match func:
        case 'delta_mean':
            try:
                ret = [argmax(a) if amax(a)-mean(a) >
                       reads_selection_threshold else False for a in preds]
            except:
                ret = [argmax(a) for a in preds]
        case 'min_max':
            try:
                ret = [argmax(a) if min([amax(a)-p for p in a if p != amax(a)]) >
                       reads_selection_threshold else False for a in preds]
            except:
                ret = [argmax(a) for a in preds]
        case 'delta_sum':
            try:
                ret = [argmax(a) if amax(a) > (sum(a)-amax(a)) + reads_selection_threshold
                       else False for a in preds]
            except:
                ret = [argmax(a) for a in preds]
        case _:
            ret = [argmax(a) for a in preds]
    if len(ret) == 0:
        raise ValueError(
            "There's no read to evaluate, your entry data might be broken.")
    elif not all(p == False for p in ret):
        my_output_msg(
            f"All reads with a threshold inferior at {reads_selection_threshold} for function {func} have been purged.")
        return ret
    else:
        return softmax_from_prediction(preds, reads_selection_threshold-0.05, func)


def softpred_from_prediction(preds, sample_name: str, clade: str, determined: str, inverted_map: dict):
    lsc = [i for i in range(len(inverted_map))] + [len(inverted_map)]
    df = pd.DataFrame(columns=lsc)
    thsh = [0.05, 0.1, 0.25, 0.5]
    # plt.style.use('seaborn-deep')

    softmax = Counter(a for a in softmax_from_prediction(
        preds, 0, 'none') if not isinstance(a, bool))

    false_reads = Counter(a for a in softmax_from_prediction(
        preds, 0, 'none') if isinstance(a, bool))[False]
    softmax[len(inverted_map)] = false_reads

    ser = pd.Series(data=softmax, index=list(softmax.keys()),
                    name=f"without filtering")

    df = df.append(ser).fillna(0)
    for t in thsh:
        for func in ['delta_mean']:  # , 'min_max', 'delta_sum'
            softmax = Counter(a for a in softmax_from_prediction(
                preds, t, func) if not isinstance(a, bool))
            false_reads = Counter(a for a in softmax_from_prediction(
                preds, t, func) if isinstance(a, bool))[False]
            softmax[len(inverted_map)] = false_reads

            ser = pd.Series(data=softmax, index=list(softmax.keys()),
                            name=f"{func}:{t}")

            df = df.append(ser).fillna(0)
    df = df.transpose()
    ax = df.plot(kind='bar', figsize=(20, 6),
                 ylabel='Reads count', rot=90, colormap='cividis')
    ax.set_xticks([i for i in range(len(inverted_map))]+[len(inverted_map)])
    ax.set_xticklabels([f"{inverted_map[str(i)]}"
                        for i in range(len(inverted_map))] + ['Rejected'])
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.13), ncol=5)
    plt.savefig(
        f"output/{sample_name}/{clade}_{determined}_softprob.svg", bbox_inches='tight')


def plot_tree_model(bst: xgb.Booster, job_name: str, classif_level: str, sp_determined: str, ksize: int) -> None:
    """Plots out one sample tree to explore params, and saves it as a .png file

    Args:
        bst (xgb.Booster): a pre-trained booster
        job_name (str): name of job, for naming purposes
        classif_level (str): level of classif we're working at
        sp_determined (str): previous level we've determined
    """
    tree = xgb.to_graphviz(bst)
    #mod_to_tree(tree, ksize, job_name, classif_level, sp_determined)


#################### CORE ####################


@my_function_timer("Building model")
def make_model(job_name: str, exclude: list[str], path: str, classif_level: str, db_name: str, sp_determined, model_parameters: dict, number_rounds: int, sample_name: str, bool_temporary):
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
        path, classif_level, 'train', db_name, sp_determined, sample_name, exclude)
    #dtrain = xgb.DMatrix(data=X, label=y)
    bst = modelisation(dtrain, model_parameters, number_rounds)
    save_model(bst, path, classif_level, db_name,
               sp_determined, job_name, bool_temporary)
    save_params(bst, path, classif_level, db_name, sp_determined)


@my_function_timer("Plotting feature importance")
def make_testing(exclude: list[str], path_to_save, size_kmer, job_name, sp_determined, path, db_name, classif_level, class_count, model_parameters, number_rounds) -> tuple:
    """_summary_

    Args:
        size_kmer (int): ksize
        job_name (str): name, for naming purposes
        sp_determined (str): previous level of classification
        path (str): path to file
        db_name (str): path to database
        classif_level (str): current classification level we're working at
        class_count (int): number of classes inside current division
        model_parameters (dict): a set of parameters
        number_rounds (int): number of boosting rounds for validation

    Returns:
        tuple: tuple containing :
            [0] -> accuracy score
            [1] -> cross-validation
            [2] -> parameters estimations
    """
    bst = xgb.Booster()
    global_scope()
    mat = load_xgboost_data(
        path, classif_level, 'train', db_name, sp_determined, job_name, exclude)
    xgb_cv = xgb.cv(dtrain=mat, params=init_parameters(class_count), nfold=3,
                    num_boost_round=20, early_stopping_rounds=10, metrics="auc", as_pandas=True, seed=123)
    bst = modelisation(mat, model_parameters, number_rounds)
    mapped = {
        f"{recode_kmer_4(str(k[1:]),size_kmer)}": v for k, v in bst.get_score(importance_type='gain').items()}
    if mapped != {}:
        plot_features(path_to_save, mapped, job_name,
                      classif_level, sp_determined)
    plot_tree_model(bst, job_name, classif_level, sp_determined, size_kmer)
    return xgb_cv
