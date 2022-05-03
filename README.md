# WISP
XGBoost implementation for determining bacteria taxas.

![img](https://github.com/Tharos-ux/wisp/blob/master/ensifer.png)

This Python program is meant to determine to which taxa a bacteria is belonginig to, from long reads (>10.000bp), solely based upon alignment-free methods. As of now, the main focus is upon kmers proportions.
Currently, five levels of taxa are implemented : **domain**, **phylum**, **group**, **order** and **family**.
This tool is requiring some reference genomes, which it will index, to create a XGBoost model.
It has been trained in its initial state over the 143 genomes assembly (see https://espace.library.uq.edu.au/view/UQ:411283)
Once a model finished at a given taxa level, it aims to do another iteration from previous results, dividing non-matching reference genomes.

# Requirements

Requires a Python 3.10 version, and all packages listed in the *requirements.txt* file.
Once your python envrionnement or installation is up, you can install quickly all dependencies with :

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

**Optional arguments :**

+ -f, -file > Specifies a file name inside output unk folder to test only this file and not the full folder

In result, you can call the file alternatively with :

```bash
python main.py db_name wisp_params.json my_sample_job -f myfnafile.fna
```

**Exploration of parameters :**

You are free to change those parameters to better fit your current classification requirements.
Code and comments are self explanatory enough, so see this snippet of code for documentation.

```python
params_job: dict = {
    # 'taxa' : [kmer_size, reads_size, subsampling_depth]
    # params for your database here
        'domain_ref': [5, 10000, 100],
        'phylum_ref': [5, 10000, 500],
        'group_ref': [5, 10000, 750],
        'order_ref': [5, 10000, 750],
        'family_ref': [5, 10000, 500],
    # params for for your sample here
        'domain_sample': [5, 10000, 500],
        'phylum_sample': [5, 10000, 750],
        'group_sample': [5, 10000, 750],
        'order_sample': [5, 10000, 750],
        'family_sample': [5, 10000, 500],
    # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/Genomes/",
    # 'output' : output for database
        'output': "data/",
    # parameters for exploration and algorithm
        'threshold': 0.1,
        'nb_boosts': 10,
    # parameters regarding results : show full takes time
        'full_test_set': True
}
```

You can modify the raw .json file, or rather update the dict contained in script *wisp_lib/parameters_init.py* and re-generates the file with :

```bash
python wisp_lib/parameters_init.py
```

# Results

Results are currently exported into a set of plots, a html file and a json file, in order to easily overview data or fetch it to another program.

See some samples here :

![img](https://github.com/Tharos-ux/wisp/blob/main/test_clostridium.png)

![img](https://github.com/Tharos-ux/wisp/blob/main/test_pseudomonas.png)

![img](https://github.com/Tharos-ux/wisp/blob/main/test_rhodococcus.png)
