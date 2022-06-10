from pdflatex import PDFLaTeX
from string import Template
from datetime import date


def teXstr(entry: str) -> str:
    """Escapes all underscore in tex string

    Args:
        entry (str): a TeX string

    Returns:
        str: a TeX string with escaped underscores
    """
    return str(entry).replace('_', '\\_')


def tex_properties(name_of_sample: str, name_of_read: str) -> str:
    """head of TeX doc (packages and title of doc)

    Args:
        name_of_sample (str): for naming purposes

    Returns:
        str: head and properties + packages used
    """
    vars_title: list[str] = [Template('\\title{JOB : $name\\\\[0.2em]\\smaller{}READ ID : $readname}').substitute(
        name=name_of_sample, readname=name_of_read), Template('\\date{$day}').substitute(day=date.today())]
    return Template('\\documentclass[12pt]{article}\n\\usepackage[a4paper, total={6in, 8in}]{geometry}\n\\usepackage[utf8]{inputenc}\n\\usepackage[clean]{svg}\n\\usepackage{graphicx}\n\\usepackage{caption}\n\\usepackage{float}\n\\usepackage{flafter}\n\\usepackage{hyperref}\n$ref').substitute(ref='\n'.join(vars_title))


def header(job_name: str, threshold: float, reads_ratio: float, number_subreads: int, version: str, nb_bp: int) -> str:
    """creates head of doc (abstract and global results)

    Args:
        job_name (str): name of task
        threshold (float): exploration percentage
        reads_ratio (float): ratio of cleared reads

    Returns:
        str: header for doc in TeX
    """
    abstract = Template('Sample has been splitted in $number distinct lectures over $nucleotids sequenced nucleotides.\\\\\nExplored hypothesis are all above $percentage percent of attributed reads.\\\\\nAll explorations have been made within a significance range of [0, $ratio[.\\\\\nThis report was produced with WISP version $version. \\\\\nYou may get source code from \\url{https://github.com/Tharos-ux/wisp}').substitute(
        percentage=int(threshold*100), ratio=reads_ratio, number=number_subreads, version=version, nucleotids=nb_bp)
    figs: str = subfigure([f"{job_name}_pie_merge.svg",
                          f"{job_name}_tree.svg"], f"Global data for {job_name}", 0.4)
    return Template('\\begin{abstract}\n\\begin{sloppypar}\n$abstract\n\\end{sloppypar}\n\\end{abstract}$figures').substitute(abstract=abstract, figures=figs)


def tex_closure() -> str:
    """Computes the end of doc

    Returns:
        str: end of TeX doc
    """
    return Template('\\end{document}').substitute()


def subfigure(list_of_figures: list[str], caption: str, w: float) -> str:
    """Creates some LaTeX subfigures with all the plots required

    Args:
        list_of_figures (list[str]): all the graphs you want included (relative paths)
        caption (str): caption for figure

    Returns:
        str: a LaTeX string with all figures
    """
    subfigure: Template = Template(
        '\\includesvg[width=$width\\textwidth]{$graph}')
    figure: Template = Template(
        '\\begin{figure}[h]\n\\centering\n$allfigures\n\\caption{$legend}\n\\label{$label}\n\\end{figure}')
    subfigures_collection: str = '\n'.join(
        [subfigure.substitute(graph=f, width=w) for f in list_of_figures])
    return figure.substitute(allfigures=subfigures_collection, legend=caption, label=f"f-{caption}")


def table(list_of_elements: list[list], caption: str) -> str:
    """Format data into a TeX table

    Args:
        list_of_elements (list[list]): list of list, [[line],[line],...]
        caption (str): a legend for the table

    Returns:
        str: a LaTeX string defining a 2D array
    """
    if not isinstance(list_of_elements, list):
        return ""
    linear_table = ' \\\\\n'.join(
        [' & '.join([str(e) for e in elt]) for elt in list_of_elements])
    return Template('\\begin{table}[htp]\n\\begin{tabular}{$size}\n$contents\n\\end{tabular}\n\\caption*{$label}\n\\end{table}').substitute(label=caption, contents=linear_table, size=''.join(['l' for _ in range(len(list_of_elements[0]))]))


def render_output(texfile: str):
    """Renders a LaTeX pdf output

    Args:
        texfile (str): path to a .tex file
    """
    pdfl = PDFLaTeX.from_texfile(f"{texfile}.tex")
    pdf, log, completed_process = pdfl.create_pdf(
        keep_pdf_file=True, keep_log_file=True)


def save_tex_file(data: str, path: str) -> None:
    """Given a .tex string, saves a .tex file

    Args:
        data (str): .tex formatted string
        path (str): output path to save file
    """
    with open(f"{path}.tex", "w") as writer:
        writer.write(data)


def make_doc(path_for_read: str, job_name: str, params: dict, taxas_levels: list[str], reports: dict, test_results: dict, test_mode: bool, threshold: float, reads_ratio: float, number_subreads: int, name_of_read: str, version, nb_bp: int) -> None:
    """Command to create and save the report

    Args:
        job_name (str): name of job, for identification and naming purposes
        params (dict): all params utilised by the algorithm
        taxas_levels (list[str]): levels inside classification
        reports (dict): all data outputted from computation
        test_results (dict): _description_
        test_mode (bool): _description_
        threshold (float): _description_
        reads_ratio (float): _description_
    """
    tex_string: str = teXstr(Template('$header\n\\begin{document}\n\\maketitle\n$head\n$docstring\n$footer').substitute(header=tex_properties(job_name, name_of_read), head=header(
        job_name, threshold, reads_ratio, number_subreads, version, nb_bp), docstring=generate_core_tex(params, taxas_levels, reports, test_results, test_mode), footer=tex_closure()))
    save_tex_file(
        tex_string, f"{path_for_read}{job_name}_{name_of_read.replace(' ','')}")
    # render_output(path) # bugged for now ; render error with /tmp


def page_break() -> str:
    """Calls for flush and newpage

    Returns:
        str: a LaTeX command for breaking to next page
    """
    return '\\clearpage\n\\pagebreak[4]'


def make_title(text: str) -> str:
    """Creates a unnumbered title for our report

    Args:
        text (str): content of the title

    Returns:
        str: a section type TeX title
    """
    return Template('\\section*{$title}').substitute(title=text)


def global_sample_report(global_path: str, job_name: str, results: dict) -> None:
    tex_string: str = teXstr(Template('$header\n\\begin{document}\n\\maketitle\n$docstring\n$footer').substitute(
        header=tex_properties(job_name, "results_all_reads"), docstring=generate_core_global(results), footer=tex_closure()))
    save_tex_file(
        tex_string, f"{global_path}{job_name}")


def generate_core_global(results):
    return table([[k, v] for k, v in results.items()], f"Individual reads results listing")


def generate_core_tex(params: dict, taxas_levels: list[str], reports: dict, test_results: dict, enlable_supplementary_estimators: bool) -> str:
    """Generate the core of the report

    Args:
        params (dict): parameters used for the booster
        taxas_levels (list[str]): list of levels used for classification
        reports (dict): results of the identification and such
        test_results (dict): results of the estimators of the booster
        enlable_supplementary_estimators (bool): tells if more graphs are needed

    Returns:
        str: a LaTeX string containing all the core of the report
    """
    # table for booster parameters
    table_params: str = table(
        [[k, v] for k, v in params.items()], 'Algorithm parameters')

    # figures & estimators for each level
    titles: list[str] = []
    levels_figures: list[str] = []
    supplementary_figures: list[str] = []
    supplementary_figures_2: list[str] = []
    levels_estimators: list[str] = []
    all_hypothesis = [['None'], reports['Possible for domain'], reports['Possible for phylum'],
                      reports['Possible for group'], reports['Possible for order'], reports['Possible for family']]
    for i, clade in enumerate(taxas_levels):
        for h in all_hypothesis[i]:
            levels_figures.append(subfigure([f"{clade}_{h}_confusion_matrix.svg", f"{clade}_{h}_graph_reads.svg"],
                                  f"Level {clade} for hypothesis {h}" if h != 'None' else f"Level {clade}", 0.8))
            supplementary_figures.append(
                subfigure([f"{clade}_{h}_boosting_results.svg", f"{clade}_{h}_feature_importance.svg"], f"Supplementary material for {h}", 0.8) if enlable_supplementary_estimators else '\n')
            supplementary_figures_2.append(
                subfigure([f"{clade}_{h}_softprob.svg", f"{clade}_{h}_proba_reads.svg"], f"Supplementary material for {h}", 0.8) if enlable_supplementary_estimators else '\n')
            table_res = [[str(k)+' : '+str(v)]
                         for k, v in test_results[f"{clade}_{h}"].items()]
            titles.append(
                f"Level {clade} for hypothesis {h}" if h != 'None' else f"Level {clade}")
            levels_estimators.append(table(table_res, f"Estimators for {h}"))
    core_results: str = '\n'.join(
        [f"{page_break()}\n{make_title(titles[i])}\n{fig}\n{levels_estimators[i]}\n{supplementary_figures[i]}\n{supplementary_figures_2[i]}" for i, fig in enumerate(levels_figures)])
    return '\n'.join([table_params, core_results])
