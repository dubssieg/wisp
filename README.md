# WISP : XGBoost implementation for bacteria families identification from long reads

![img](https://github.com/Tharos-ux/wisp/blob/master/preview/Ente_malo.png)

This Python program is meant to determine to which taxa a bacteria is belonginig to, from long reads (>10.000bp), solely based upon alignment-free methods. As of now, the main focus is upon kmers proportions.
Currently, five levels of taxa are implemented : **domain**, **phylum**, **group**, **order** and **family**.
This tool is requiring some reference genomes, which it will index, to create a XGBoost model.
It has been trained in its initial state over the 143 genomes assembly (see https://espace.library.uq.edu.au/view/UQ:411283)
Once a model finished at a given taxa level, it aims to do another iteration from previous results, excluding non-matching reference genomes.
The core functionnalities relies on a class probabiliy attribution to discriminate reads that might not be good indicators for our specie to be determined. As many other options, you can choose the ratio and the selection function to suit best your biological context.

# Requirements

Requires a Python 3.10 version, and all packages listed in the *requirements.txt* file.
Once your python envrionnement or installation is up, you can install quickly all dependencies with :

```bash
pip install -r build/requirements.txt
```

Also requires a Graphviz installation in order to perform trees rendering.

# Launching WISP with wisp.py

This way is best suited for huge amounts of genomes computations : it will iterate through your unknown genomes and create separate reports for each one.
Will be later on modified to a bash script in order to add this tool more easily into pipelines.

```bash
python wisp.py [-t number]
```

**Optional arguments :**

+ -t, --multithreading > you can imput an int to define the number of parallel processes you want to spawn (recommanded : your number of cores)

# Launching WISP with main.py

This way is more suited for single genome analysis, or intensive plotting for one result.
From command line, main class can be used with the following arguments :

**Required arguments :**

```bash
python main.py db_name wisp_params.json my_sample_job
```

+ "db_name" > Name of database, will buid it if ont exists (type=str)
+ "params" > Path to a .json parameters file. See example.
+ "job_name" > Name of job (type=str)

**Optional arguments :**

+ -f, -file > Specifies a file name inside output unk folder to test only this file and not the full folder. If no file is specified, all files will be merged as one single sample.

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
        'domain_ref': [4, 10000, 50],
        'phylum_ref': [4, 10000, 250],
        'group_ref': [4, 10000, 375],
        'order_ref': [4, 10000, 375],
        'family_ref': [4, 10000, 250],
        # params for your sample here
        'domain_sample': [4, 10000, 500],
        'phylum_sample': [4, 10000, 750],
        'group_sample': [4, 10000, 750],
        'order_sample': [4, 10000, 750],
        'family_sample': [4, 10000, 500],
        # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/Genomes/",
        # 'output' : output for database
        'output': "data/",
        # parameters for exploration and algorithm
        'threshold': 0.2,
        'nb_boosts': 12,
        # parameters regarding results
        'full_test_set': False,
        # parameter for read selection, signifiance for softprob
        'reads_th': 0.25,
        # force rebuilding full model, when you re-use a database but you changed model parameters
        'force_model_rebuild': True
    }
```

You can modify the raw .json file, or update the dict contained in script *wisp_lib/parameters_init.py* and re-generate the file with :

```bash
python wisp_lib/parameters_init.py
```

# Force database building

One can enforce full database building (all levels and all branches) to generate full database and models for future usage. This folder then can be safely copied to another installation of WISP and be used from here.
It aims to allow people to index huge databases on reliable machines, then loading it on another less powerfull computer.

```bash
python force_build.py db_name params.json 
```

+ "db_name" > Name of database, will buid it if ont exists (type=str)
+ "params" > Path to a .json parameters file. See example.

# Results

Results are currently exported into a set of plots, a html file and a json file, in order to easily overview data or fetch it to another program.

**How to interpret a report**

For each classification level, you get a set of infos. You may get even more with `'full_test_set': True` but compute time will be higher. The first row tells the global classification result, based upon the reads selection method you choose.
At the end of the report, you may find the full parameters set you used on this particular sample, for archiving purposes.
The tree displays information on which ways were investigated, matching the threshold you defined. Nodes are classes and edges are percentages of attributed reads on those classes from previous node.

![tree_example](https://github.com/Tharos-ux/wisp/blob/master/preview/tree_example.png)

Then, you get the information by level : in base mode, you have one confusion matrix, displaying the results from a test set against the trained model ; numbers are row percentages, showing for each actual class the percentage of correctly (or not) predicted reads. Note that this info is calculated upon your selection parameters, so you can investigate which parameters works better for your database, but you should look out for overfitting!

![confusion_matrix](https://github.com/Tharos-ux/wisp/blob/master/preview/phylum_confusion_matrix.png)

The next graph tells you about reads attribution : once calculations are done, reads are attributed to a certain class.
The dashed line displays the investigation limit you set (reads percentage to consider a class)

![graph_reads](https://github.com/Tharos-ux/wisp/blob/master/preview/phylum_bacteria_graph_reads.png)

Next up, if you choose advanced plotting options, you get an additionnal set of graphs.
Firstly, you have 'Mean and standard deviation' which tells you more about the fidelity of successive boostings.

![mean_deviation](https://github.com/Tharos-ux/wisp/blob/master/preview/family_lactobacillales_boosting_results.png)

Secondly, you get a plotting of the 15 most important features for the split, to get an extra bit more info upon the signature concept which this algorithm relies on.

![features_example](https://github.com/Tharos-ux/wisp/blob/master/preview/order_bacilli_feature_importance.png)

Lastly, you get additional plots that are not included inside the .html file but are freely consultable from the output folder :

+ Reads porbability repartition
+ Reads selection across functions

# Known errors

+ Rare bug where re-running a job will not be procedeed if report is opened in webbrowser and folder was deleted few seconds ago
+ Rare bug in wisp multithreading mode where the sample file would not be created