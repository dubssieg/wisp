# WISP
XGBoost implementation for determining bacteria taxas.

![img](https://github.com/Tharos-ux/wisp/blob/main/test_lactiplantibacillus.png)

This Python program is meant to determine to which taxa a bacteria is belonginig to, from long reads (>10.000bp), solely based upon alignment-free methods. As of now, the main focus is upon kmers proportions.
Currently, four levels of taxa are implemented : **domain**, **phylum**, **order** and **family**.
This tool is requiring some reference genomes, which it will index, to create a XGBoost model.
It has been trained in its initial state over the 143 genomes assembly (see https://espace.library.uq.edu.au/view/UQ:411283)
Once a model finished at a given taxa level, it aims to do another iteration from previous results, dividing non-matching reference genomes.

# Requirements

Install all dependencies with :

```bash
pip install -r requirements.txt
```

# Parameters

From command line, main class can be used with the following arguments :

**Required arguments :**

```bash
python main.py db_name wisp_params.json my_sample_job
```

+ "db_name" > Name of database, will buid it if ont exists (type=str)
+ "params" > Path to a .json parameters file. See example.
+ "job_name" > Name of job (type=str)

# Results

Results are currently exported into a set of plots, a html file and a json file, in order to easily overview data or fetch it to another program.

See some samples here :

![img](https://github.com/Tharos-ux/wisp/blob/main/test_clostridium.png)

![img](https://github.com/Tharos-ux/wisp/blob/main/test_pseudomonas.png)

![img](https://github.com/Tharos-ux/wisp/blob/main/test_rhodococcus.png)
