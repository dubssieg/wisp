from collections import Counter
from multiprocessing.dummy import Array
from numpy import array, diag
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sns
from pandas import DataFrame, Series, crosstab


def reads_species_plotter(predicitions, sample_name: str, inverted_map: dict, clade: str, determined: str, threshold: float) -> None:
    """Plots out reads repartition

    Args:
        predicitions (xgboost preds): full dataset of predictions
        sample_name (str): name of job, used for naming purposes
        inverted_map (dict): mapping from classes to clades
        clade (str): level of classification we're working at
        determined (str): previous hypothesis of determination
        threshold (float): to plot a line, threshold for exploration options
    """
    datas = dict(Counter(predicitions))
    keys_sorted = sorted([k for k in datas.keys()])
    x = array(keys_sorted)
    y = array([datas[k] for k in keys_sorted])
    maxm = sum(y)
    figure = plt.figure(figsize=(9, 6))
    graph = figure.add_axes([0.1, 0.2, 0.8, 0.6])
    bars = graph.bar(x, y, color='black')
    plt.title(
        f"Reads accross {clade} for {sample_name}")
    plt.ylabel("Read counts")
    plt.xticks(rotation=70, fontsize=6, ha='right')
    graph.set_xticks([i for i in range(len(inverted_map))])
    graph.set_xticklabels([f"{inverted_map[str(i)]}"
                          for i in range(len(inverted_map))])
    graph.bar_label(bars)
    plt.axhline(y=threshold*maxm, xmin=0, color="#465065",
                linestyle='dashed', linewidth=1)
    plt.savefig(f"output/{sample_name}/{clade}_{determined}_graph_reads.png")


def plot_all_reads(matrix, sample_name: str, inverted_map: dict, clade: str, determined: str):
    figure = plt.figure(figsize=(10, 6))
    graph = figure.add_axes([0.1, 0.2, 0.8, 0.6])
    for serie in matrix:
        plt.plot(serie, color="#465065")
    plt.title(
        f"Probabilities accross {clade} for {sample_name}")
    plt.ylabel("Probabilities")
    plt.xticks(rotation=70)
    graph.set_xticks([i for i in range(len(inverted_map))])
    graph.set_xticklabels([f"{inverted_map[str(i)]}"
                          for i in range(len(inverted_map))])
    plt.savefig(f"output/{sample_name}/{clade}_{determined}_proba_reads.png")


def compare_test(test_classes, test_preds, inverted_map: dict, sample_name: str, clade: str, determined: str):
    """Calling for estimators and plotting of confusion matrix

    Args:
        test_classes (array): true classes we're dealing wih for our test set
        test_preds (array): predicted classes for test set
        inverted_map (dict): mapping between ints and class labels
        sample_name (str): used for naming purposes
        clade (str): level we're working at
        determined (str): upper level already determined

    Returns:
        _type_: _description_
    """
    test_preds = [t for t in test_preds]
    cm = pandas_confusion(test_classes, test_preds,
                          inverted_map)
    plot_pandas(cm, sample_name, clade, determined)
    return text_classification_report(test_classes, test_preds, inverted_map)


def text_classification_report(test_classes, test_preds, inverted_map: dict) -> dict:
    """Gives out a report for all test conditions with estimators

    Args:
        test_classes (array): true classes we're dealing wih for our test set
        test_preds (array): predicted classes for test set
        inverted_map (dict): mapping between ints and class labels

    Returns:
        dict: _description_
    """
    cls = classification_report(test_classes, test_preds, output_dict=True)
    for key, value in list(cls.items()):
        if key in inverted_map.keys():
            cls[inverted_map[key]] = cls.pop(key)
        if isinstance(value, dict):
            for k, v in value.items():
                if isinstance(v, float):
                    value[k] = round(value[k], 2)
        if isinstance(value, float):
            value = round(value, 2)
    return cls


def pandas_confusion(test_classes, test_preds, inverted_map: dict) -> DataFrame:
    """Creates the dataframe used to plot a confusion matrix from test datas

    Args:
        test_classes (array): true classes we're dealing wih for our test set
        test_preds (array): predicted classes for test set
        inverted_map (dict): mapping between ints and class labels

    Returns:
        pd.DataFrame: confusion matrix
    """
    test_classes = [inverted_map[str(int(t))] for i, t in enumerate(
        test_classes) if not isinstance(test_preds[i], bool)]
    test_preds = [inverted_map[str(int(t))]
                  for t in test_preds if not isinstance(t, bool)]
    data = {'y_Actual': test_classes, 'y_Predicted': test_preds}
    df = DataFrame(data, columns=['y_Actual', 'y_Predicted'])
    return crosstab(df['y_Actual'], df['y_Predicted'], rownames=[
        'Actual'], colnames=['Predicted'])


def plot_pandas(cm: DataFrame, sample_name: str, clade: str, determined: str, cmap: str = 'bone') -> dict:
    """Plots the confusion matrix for test data at given level

    Args:
        cm (pd.DataFrame): a pandas df that contains a confusion matrix
        sample_name (str): used for naming purposes
        clade (str): level we're working at
        determined (str): upper level already determined
        cmap (str, optional): set of colors for the heatmap. Defaults to 'bone'.

    Returns:
        dict : number of true and false classifications for clade
    """
    plt.figure(figsize=(7, 6))
    ax = plt.axes()
    try:
        diag_axis = Series(diag(cm)).sum()
        full_set = cm.to_numpy().sum()
        ax.set_title(
            f"Confusion matrix for {clade} level. R={round(diag_axis/(full_set-diag_axis),2)}")
    except:
        ax.set_title(
            f"Confusion matrix for {clade} level.")
    # percentage
    cm = cm.div(cm.sum(axis=1), axis=0) * 100
    sns.heatmap(cm, annot=True, cmap=cmap, fmt='.0f', linewidths=0.5, ax=ax)
    plt.yticks(fontsize=5)
    plt.xticks(fontsize=5)
    plt.savefig(
        f"output/{sample_name}/{clade}_{determined}_confusion_matrix.png")


def plot_boosting(df, sample_name, clade, determined, number_rounds):

    df.reset_index()
    X = df.index
    mean_train = df.iloc[:, 0]
    mean_test = df.iloc[:, 2]
    sd_train = df.iloc[:, 1]
    sd_test = df.iloc[:, 3]

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
    fig.suptitle(
        f"Mean and standard deviation for {clade} with hypothesis {determined}")
    axs[0].plot(X, mean_train, color='#53668f', label='train')
    axs[0].plot(X, mean_test, color='#07142f', label='test')
    axs[0].vlines(number_rounds, ymin=0.9, ymax=1, color="#465065",
                  linestyle='dashed', linewidth=1)
    axs[0].legend(loc='lower right')
    axs[1].plot(X, sd_train, color='#53668f', label='train')
    axs[1].plot(X, sd_test, color='#07142f', label='test')
    axs[1].vlines(number_rounds, ymin=0, ymax=0.01, color="#465065",
                  linestyle='dashed', linewidth=1)
    axs[1].legend(loc='upper right')
    plt.savefig(
        f"output/{sample_name}/{clade}_{determined}_boosting_results.png")


def plot_features(datas, job_name, classif_level, sp_determined):

    keys = list(datas.keys())
    values = list(datas.values())

    data = DataFrame(data=values, index=keys, columns=[
        "score"]).sort_values(by="score", ascending=False)
    data.astype(float).nlargest(15, columns="score").plot(
        kind='barh', color='#465065', figsize=(7, 6))
    plt.savefig(
        f"output/{job_name}/{classif_level}_{sp_determined}_feature_importance.png")
