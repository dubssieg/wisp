import pygraphviz as PG


def tree_render(results: dict, job_name: str, path: list) -> None:
    """Renders a classification tree with pygraphviz engine

    Args:
        results (dict): dictionnary that contains all the results of all the unk runs
        job_name (str): name of the job, to define the output
        path (list): list of clades we're assuming is the correct one
    """
    root = ["None (-)"]
    tree = PG.AGraph(directed=False, strict=True)
    unpacking(tree, root, results, path)
    tree.layout(prog='dot')
    tree.draw(f"output/{job_name}/{job_name}_tree.png")


def unpacking(tree: PG.AGraph, root: list, datas: dict, path: list) -> None:
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
                    tree.add_node(f"{rt}", color="#a7b7d9")
            unpacking(tree, new_root, datas, path)
        except:
            pass
