from collections import Counter
from numpy import array, diag
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
import seaborn as sns
from pandas import DataFrame, Series, crosstab


def reads_species_plotter(path_to_save: str, predicitions, sample_name: str, inverted_map: dict, clade: str, determined: str, threshold: float, sampling_number) -> None:
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
    rejected = sampling_number - sum([datas[k] for k in keys_sorted])
    x = array(keys_sorted + [max([int(i) for i in inverted_map.keys()])+1])
    templist = [datas[k] for k in keys_sorted] + [rejected]
    y = array(templist)
    figure = plt.figure(figsize=(9, 6))
    graph = figure.add_axes([0.1, 0.2, 0.8, 0.6])
    colors = ['#454261' if datas[k] != max(
        [datas[k] for k in keys_sorted]) else '#1c1a2c' for k in keys_sorted] + ['#594966']
    bars = graph.bar(x, y, color=colors)
    plt.title(
        f"Reads accross {clade} for {sample_name}")
    plt.ylabel("Read counts")
    plt.xticks(rotation=70, fontsize=9, ha='right')
    tl = [i for i in range(len(inverted_map))
          ]+[len(inverted_map)]
    graph.set_xticks(tl)
    tk = [f"{inverted_map[str(i)]}"
          for i in range(len(inverted_map))] + ["Rejected"]
    graph.set_xticklabels(tk)
    graph.bar_label(bars)
    plt.axhline(y=threshold*sampling_number, xmin=0, color="#465065",
                linestyle='dashed', linewidth=1)
    plt.savefig(
        f"{path_to_save}{clade}_{determined}_graph_reads.png", bbox_inches='tight')


def plot_all_reads(matrix, sample_name: str, inverted_map: dict, clade: str, determined: str) -> None:
    """Does the plotting of reads repartition across clades

    Args:
        matrix (numpy array): set of predictions softmaxed we want to plot
        sample_name (str): for naming purposes
        inverted_map (dict): mapping between numbers and clades names
        clade (str): level we're working at
        determined (str): level previously determined
    """
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
    plt.savefig(
        f"output/{sample_name}/{clade}_{determined}_proba_reads.png", bbox_inches='tight')


def plot_pie_merge(path_to_save: str, prediction_reads: dict[str, float], sample_name: str) -> None:
    values, labels = [], []
    for k, v in prediction_reads.items():
        values.append(v)
        labels.append(k)
    if sum(values) < 1:
        values.append(1-sum(values))
        labels.append('Rejected')
    figure = plt.figure(figsize=(5, 4))
    plt.pie(values, labels=labels)  # , rotatelabels=True
    plt.axis('equal')

    plt.title(
        f"Raw reads predicitions for {sample_name}")
    plt.savefig(
        f"{path_to_save}{sample_name}_pie_merge.png", bbox_inches='tight')


def compare_test(path_to_save: str, test_classes, test_preds, inverted_map: dict, sample_name: str, clade: str, determined: str) -> dict:
    """Calling for estimators and plotting of confusion matrix

    Args:
        test_classes (array): true classes we're dealing wih for our test set
        test_preds (array): predicted classes for test set
        inverted_map (dict): mapping between ints and class labels
        sample_name (str): used for naming purposes
        clade (str): level we're working at
        determined (str): upper level already determined

    Returns:
        dict: report with estimators
    """
    test_preds = [t for t in test_preds]
    cm = pandas_confusion(test_classes, test_preds,
                          inverted_map)
    plot_pandas(path_to_save, cm, sample_name, clade, determined)
    return text_classification_report(test_classes, test_preds, inverted_map)


def text_classification_report(test_classes, test_preds, inverted_map: dict) -> dict:
    """Gives out a report for all test conditions with estimators

    Args:
        test_classes (array): true classes we're dealing wih for our test set
        test_preds (array): predicted classes for test set
        inverted_map (dict): mapping between ints and class labels

    Returns:
        dict: report with estimators
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


def plot_pandas(path_to_save: str, cm: DataFrame, sample_name: str, clade: str, determined: str, cmap: str = 'bone') -> None:
    """Plots the confusion matrix for test data at given level

    Args:
        cm (pd.DataFrame): a pandas df that contains a confusion matrix
        sample_name (str): used for naming purposes
        clade (str): level we're working at
        determined (str): upper level already determined
        cmap (str, optional): set of colors for the heatmap. Defaults to 'bone'.
    """
    number_digits: int = 2
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
    sns.heatmap(cm, annot=True, cmap=cmap,
                fmt=f".{number_digits}f", linewidths=0.5, ax=ax)
    # plt.yticks(fontsize=9)
    # plt.xticks(fontsize=9)
    plt.savefig(
        f"{path_to_save}{clade}_{determined}_confusion_matrix.png", bbox_inches='tight')


def plot_boosting(df: DataFrame, sample_name: str, clade: str, determined: str, number_rounds: int) -> None:
    """Plots the gain from multiple boostings

    Args:
        df (DataFrame): datas of each series
        sample_name (str): name of job, naming purposes
        clade (str): level we're working at
        determined (str): previous level of classification
        number_rounds (int): number of successive boostings that were done
    """
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
        f"output/{sample_name}/{clade}_{determined}_boosting_results.png", bbox_inches='tight')


def plot_features(path_to_save: str, datas, job_name: str, classif_level: str, sp_determined: str) -> None:
    """Plots the repartition of features in gain mode which were used in our model

    Args:
        datas (Counter): a count of features
        job_name (str): for naming purposes
        classif_level (str): level we're working at
        sp_determined (str): previous level estimation
    """
    keys: list = list(datas.keys())
    values: list = list(datas.values())
    nb_features: int = 16

    data = DataFrame(data=values, index=keys, columns=[
        "score"]).sort_values(by="score", ascending=False)
    data.astype(float).nlargest(nb_features, columns="score").plot(
        kind='barh', color='#465065', figsize=(7, 6))
    plt.xlabel("Gain score")
    if sp_determined != 'None':
        plt.title(
            f"Top {nb_features} features for {classif_level} (knowing {sp_determined})")
    else:
        plt.title(
            f"Top {nb_features} features for {classif_level}")
    plt.savefig(
        f"{path_to_save}{classif_level}_{sp_determined}_feature_importance.png", bbox_inches='tight')
