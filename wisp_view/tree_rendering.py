from pygraphviz import AGraph


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


def branch_score(from_node: str, tree_stats: dict) -> int:
    """Returns the score of a branch, solely looking upon number of co-branches

    Args:
        from_node (str): node we're starting from
        tree_stats (dict): a dict outputted by tree_stats

    Returns:
        int: a score, 1 means perfect path, more means less than perfect
    """
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


def tree_stats(tree: AGraph) -> dict:
    """Fetches the scores for each node

    Args:
        tree (AGraph): tree we're investigating

    Returns:
        dict: a dict of scores, ranging [1,inf]
    """
    return {node: node_childrens(tree, node) for node in names_children(tree, "None (-)")}


def final_node_score(tree: AGraph, lon: list) -> dict:
    """Evaluates the score for each terminal node

    Args:
        tree (AGraph): tree we're investigating
        lon (list): a list of leaf nodes

    Returns:
        dict: a dict of score of each leaf
    """
    my_stats = tree_stats(tree)
    return {node: sum([my_stats[i] for i in names_predecessors(tree, node) if i != 'None (-)']) for node in lon}


def tree_evaluator(tree: AGraph, path: list[str]) -> str:
    """Main routine : calls calculations upon a tree

    Args:
        tree (AGraph): tree we're investigating into
        path (list[str]): the prediction path to check against

    Returns:
        str: a sentence that resumes the parcimony of the tree
    """
    scores = final_node_score(
        tree, [p for p in names_children(tree, 'None (-)') if p[-2] == 'f'])
    if min(scores, key=scores.get) in path:
        return f"Default path is more or equally parcimonious than any other."
    else:
        return f"Though is it not the final guess, path leading to {min(scores, key=scores.get)[:-4]} is the most parcimonious."


def tree_render(results: dict, job_name: str, path: list) -> str:
    """Renders a classification tree with pygraphviz engine

    Args:
        results (dict): dictionnary that contains all the results of all the unk runs
        job_name (str): name of the job, to define the output
        path (list): list of clades we're assuming is the correct one
    """
    root = ["None (-)"]
    tree = AGraph(directed=False, strict=True)
    unpacking(tree, root, results, path)
    tree.layout(prog='dot')
    tree.draw(f"output/{job_name}/{job_name}_tree.png")
    return tree_evaluator(tree, path)


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
