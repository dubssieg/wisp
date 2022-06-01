from re import template
from pdflatex import PDFLaTeX
from string import Template
from datetime import date


def tex_properties(name_of_sample: str) -> str:
    vars_title: list[str] = [Template('\\title{$name}').substitute(
        name=name_of_sample), Template('\\date{$day}').substitute(day=date.today())]
    return Template('\\documentclass[12pt]{article}\n\\usepackage[a4paper, total={6in, 8in}]{geometry}\n\\usepackage[utf8]{inputenc}\n\\usepackage{graphicx}\n\\usepackage{hyperref}\n$ref').substitute(ref='\n'.join(vars_title))


def subfigure(list_of_figures: list[str], caption: str) -> str:
    """Creates some LaTeX subfigures with all the plots required

    Args:
        list_of_figures (list[str]): all the graphs you want included (relative paths)
        caption (str): caption for figure

    Returns:
        str: a LaTeX string with all figures
    """
    subfigure: Template = Template(
        '\\includegraphics[width=0.40\\textwidth]{$graph}')
    figure: Template = Template(
        '\\begin{figure}[h]\n\\centering\n$allfigures\n\\caption{$legend}\n\\label{fig:mesh1}\n\\end{figure}')
    subfigures_collection: str = '\n'.join(
        [subfigure.substitute(graph=f) for f in list_of_figures])
    return figure.substitute(allfigures=subfigures_collection, legend=caption)


def table(list_of_elements: list[list], caption: str) -> str:
    linear_table = ' \\\\\n'.join(
        [' & '.join([str(e) for e in elt]) for elt in list_of_elements])
    return Template('\\begin{table}[]\n\\begin{tabular}{$size}\n$contents\n\\end{tabular}\n\\end{table}').substitute(contents=linear_table, size=''.join(['l' for _ in range(len(list_of_elements[0]))]))


def render_output(texfile: str):
    """Renders a LaTeX pdf output

    Args:
        texfile (str): path to a .tex file
    """
    pdfl = PDFLaTeX.from_texfile('my_file.tex')
    pdf, log, completed_process = pdfl.create_pdf(
        keep_pdf_file=True, keep_log_file=True)


def section(section_label: str, datas: str) -> str:
    return Template('\\section{$nom}\n$datas').substitute(nom=section_label, contents=datas)


def make_doc(job_name):
    return Template('$header\n\\begin{document}\n\\maketitle\n$docstring\n\\end{document}').substitute(header=tex_properties((job_name)), docstring='')


print(subfigure(['boost.png', 'ty.png', 'mlep.png'], 'Cute plot uwu'))
print(make_doc('my_sample'))
print(table([[0, 2, 3, 9, 12], [8, 3, 6, 7, 10]], ''))
