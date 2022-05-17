"From a dict, generates HTML report of results"
#import webbrowser
from os import listdir


def cute_value(val):
    if isinstance(val, dict):
        return '<br>'.join([f"{key} : {value}" for key, value in val.items()])
    else:
        return val


def gen_html_report(params, job_name, dataset, reports, list_clades, test_results, threshold: float, test_status: bool, reads_ratio=None):
    doctype: str = "<!DOCTYPE html>"
    dataset = [data for data in dataset]

    style_sheet: str = """
        <style type="text/css">
        body {
        margin-top: 1.0em;
        background-color: #ffffff;
        font-family: Helvetica, Arial, FreeSans, san-serif;
        color: #000000;
        }
        h1 { font-size: 1.8em; color: #000000; }
        h2 { font-size: 1.6em; color: #000000; }
        h3 { font-size: 1.4em; color: #000000; }
        a { color: #1e36ce; text-decoration:none; }
        tt { font-family: "Courier New", Courier; }
        /* pre style removed because it will interfer with pygments */
        p { text-indent: 0px; }
        hr { border: 0; width: 80%; border-bottom: 1px solid #aaa}
        p.caption { width: 80%; font-style: normal; text-align: left; }
        hr.figure { border: 0; width: 80%; border-bottom: 1px solid #aaa}
        .alert-text-small   { font-size: 80%;  }
        .alert-text-large   { font-size: 130%; }
        .alert-text-normal  { font-size: 90%;  }
        .alert {
        padding:8px 35px 8px 14px; margin-bottom:18px;
        text-shadow:0 1px 0 rgba(255,255,255,0.5);
        border:1px solid #bababa;
        border-radius: 4px;
        -webkit-border-radius: 4px;
        -moz-border-radius: 4px;
        color: #555;
        background-color: #f8f8f8;
        background-position: 10px 5px;
        background-repeat: no-repeat;
        background-size: 38px;
        padding-left: 55px;
        width: 75%;
        }
        .row { display: flex; }
        .column { flex: 33%; }
        .center { display: block; margin-left: auto; margin-right: auto; width: 80%; }
        .container { height: 600px; margin-left: auto; margin-right: auto; margin-top: auto; margin-bottom: auto; }
        .container img { max-height: 100%; max-width: 100%; }
        div { text-align: justify; text-justify: inter-word; }
        </style>
    """

    #abstract: str = '<br>'.join([f"{k} > {v}" for k, v in reports.items()])

    def balise(name: str, content: str) -> str:
        return f"<{name}>{content}</{name}>"

    table_params_booster = f"""<div class="column"><h1 style="text-align:center">Parameters</h1></div>
        <div class="row">
            <div class="column">
                {''.join(["<p style='text-align:center'>"+k+"</p>" for k in params.keys()])}
            </div>
            <div class="column">
                {''.join(["<p style='text-align:center'>"+str(v)+"</p>" for v in params.values()])}
            </div>
        </div>
        """

    reads_ratio_text = f"All explorations have been made within a significance range of [0, {reads_ratio}[." if reads_ratio != None else ""
    # defining figures
    string = f"""
        <div class="row">
            <div class="column">
                <h2 style="text-align:center">{reports['domain']}</h2><p style="text-align:center">(domain)</p>
            </div>
            <div class="column">
                <h2 style="text-align:center">{reports['phylum']}</h2><p style="text-align:center">(phylum)</p>
            </div>
            <div class="column">
                <h2 style="text-align:center">{reports['group']}</h2><p style="text-align:center">(group)</p>
            </div>
            <div class="column">
                <h2 style="text-align:center">{reports['order']}</h2><p style="text-align:center">(order)</p>
            </div>
            <div class="column">
                <h2 style="text-align:center">{reports['family']}</h2><p style="text-align:center">(family)</p>
            </div>
        </div>
        <h2 style="text-align:center">Explored hypothesis (above {int(threshold*100)}% of reads)</h2>
        <p style="text-align:center">{reads_ratio_text}</p>
        <div class="row">
            <div class="container">
                <img src="{job_name}_tree.png" class="center">
            </div>
        </div>
        <p style="text-align:center">{reports["parcimonious_path"]}</p>
    """
    all_hypothesis = [['None'], reports['Possible for domain'], reports['Possible for phylum'],
                      reports['Possible for group'], reports['Possible for order'], reports['Possible for family']]
    for i in range(len(list_clades)):
        #list_graphs = [file for file in listdir(f"output/{job_name}/") if list_clades[i] in file]
        #set_hypotesis = set([graph.split('_')[1] for graph in list_graphs])
        # print(set_hypotesis)
        set_hypothesis = all_hypothesis[i]

        for hypotesis in set_hypothesis:

            table_res = [
                f"<div class=\"column\"><h3 style=\"text-align:center\">{k}</h3><p style=\"text-align:center\">{cute_value(v)}</p></div>" for k, v in test_results[f"{list_clades[i]}_{hypotesis}"].items()]

            temp = "<div class=\"row\">"
            for x in range(len(table_res)):
                if x % 6 != 0:
                    temp = f"{temp}{table_res[x]}"
                else:
                    temp = f"{temp}</div><div class=\"row\">"
            res = f"{temp}</div>"

            # will add supplementary plots only if required
            if test_status:
                insert_tests = f"""
                <div class="column">
                    <img src="{list_clades[i]}_{hypotesis}_boosting_results.png">
                </div>
                <div class="column">
                    <img src="{list_clades[i]}_{hypotesis}_feature_importance.png">
                </div>
                """
            else:
                insert_tests = ""

            title = f"Level {list_clades[i]} for hypothesis {hypotesis}" if hypotesis != 'None' else f"Level {list_clades[i]}"
            string = f"""{string}
            <h2 style="text-align:center">{title}</h2>
            <div class="row">
                <div class="column">
                    <img src="{list_clades[i]}_{hypotesis}_confusion_matrix.png">
                </div>
                <div class="column">
                    <img src="{list_clades[i]}_{hypotesis}_graph_reads.png">
                </div>
                {insert_tests}
            </div>
            {res}
            """

    html_string = balise("html", balise("head", balise("title", f"Reports for job {job_name}")+style_sheet)+balise(
        "body", balise("h1 style=\"text-align:center\"", "Global reports")+string+table_params_booster))

    with open(f"output/{job_name}/{job_name}_report.html", "w") as writer:
        writer.write(html_string)

    # if you want it to open automatically
    #webbrowser.open('file://' + path.realpath(f"output/{job_name}/{job_name}_report.html"))
