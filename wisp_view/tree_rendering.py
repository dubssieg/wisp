from pygraphviz import AGraph
from numpy import argmin


def node_childrens(tree: AGraph, name_of_node: str) -> int:
    """Counts children nodes of selected node

    Args:
        tree (PG.AGraph): a tree to seek node in
        name_of_node (str): a node name in tree

    Returns:
        int: number of children of given node
    """
    try:
        res = len(tree.successors(tree.get_node(name_of_node)))
        if res != 0:
            for subnode in tree.successors(tree.get_node(name_of_node)):
                res += node_childrens(tree, subnode)
        return res
    except:
        return 0


def names_children(tree: AGraph, name_of_node: str) -> list:
    """Counts children nodes of selected node

    Args:
        tree (PG.AGraph): a tree to seek node in
        name_of_node (str): a node name in tree

    Returns:
        int: number of children of given node
    """
    try:
        res = tree.successors(tree.get_node(name_of_node))
        if res != 0:
            for subnode in tree.successors(tree.get_node(name_of_node)):
                res += names_children(tree, subnode)
        return res
    except:
        return []


def names_predecessors(tree: AGraph, name_of_node: str) -> list:
    """Counts children nodes of selected node

    Args:
        tree (PG.AGraph): a tree to seek node in
        name_of_node (str): a node name in tree

    Returns:
        int: number of children of given node
    """
    try:
        res = tree.predecessors(tree.get_node(name_of_node))
        if res != 0:
            for subnode in tree.predecessors(tree.get_node(name_of_node)):
                res += names_children(tree, subnode)
        return res
    except:
        return []


def branch_score(from_node: str, tree_stats: dict):
    taxa_level = from_node[-2]
    taxa_decrease = {
        '-': 5,
        'd': 4,
        'p': 3,
        'g': 2,
        'o': 1,
        'f': 0
    }
    return tree_stats[from_node] - taxa_decrease[taxa_level]


def tree_stats(tree):
    return {node: node_childrens(tree, node) for node in names_children(tree, "None (-)")}


def final_node_score(tree: AGraph, lon: list):
    my_stats = tree_stats(tree)
    return {node: sum([my_stats[i] for i in names_predecessors(tree, node) if i != 'None (-)']) for node in lon}


def tree_evaluator(tree, path):
    my_stats = tree_stats(tree)
    for i, level in enumerate(['d', 'p', 'g', 'o']):
        fixed_list = list(my_stats.keys())
        score_to_compare, min_all_other_scores = branch_score(path[i], my_stats), [branch_score(
            item, my_stats) for item in [key for key in fixed_list if key[-2] == level]]
        if score_to_compare >= min(min_all_other_scores):
            print(
                f"At level ({level}), default path is more or equally parcimonious than any other")
        else:
            print(
                f"At level ({level}), default path is less parcimonious than {fixed_list[argmin(min(min_all_other_scores))]}")
    print(final_node_score(
        tree, [p for p in names_children(tree, 'None (-)') if p[-2] == 'f']))


def tree_render(results: dict, job_name: str, path: list) -> None:
    """Renders a classification tree with pygraphviz engine

    Args:
        results (dict): dictionnary that contains all the results of all the unk runs
        job_name (str): name of the job, to define the output
        path (list): list of clades we're assuming is the correct one
    """
    root = ["None (-)"]
    tree = AGraph(directed=False, strict=True)
    unpacking(tree, root, results, path)
    tree_evaluator(tree, path)
    tree.layout(prog='dot')
    tree.draw(f"output/{job_name}/{job_name}_tree.png")


def unpacking(tree: AGraph, root: list, datas: dict, path: list) -> None:
    """Recursive function that builds the nodes and edges of classification tree

    Args:
        tree (PG.AGraph): tree we will build in
        root (list): list of root nodes (n-1 level)
        datas (dict): dictionnary that contains all the data needed to build and label
        path (list): default selected path, output of full classifier
    """
    for elt in root:
        try:
            new_root = datas[f"Tree {elt}"]
            for rt in new_root:
                if elt not in path or rt not in path:
                    tree.add_edge(
                        f"{elt}", f"{rt}", label=f" {datas[rt]}", color="#a7b7d9")
                else:
                    tree.add_edge(
                        f"{elt}", f"{rt}", label=f" {datas[rt]}", penwidth=2)
                if rt not in path:
                    tree.add_node(f"{rt}", color="#a7b7d9", shape="box")
            unpacking(tree, new_root, datas, path)
        except:
            pass
