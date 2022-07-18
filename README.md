# WISP : Bacterial families identification from long reads, machine learning with XGBoost

[![](https://img.shields.io/badge/Build-stable-green.svg)](https://github.com/Tharos-ux/wisp/)
[![](https://img.shields.io/badge/Wiki-unfinished-red.svg)](https://github.com/Tharos-ux/wisp/wiki)
[![](https://img.shields.io/badge/Paper-preview-blue.svg)](preview/Bacterial_families_identification_from_long_reads___machine_learning_with_XGBoost_.pdf)
[![](https://img.shields.io/badge/Licence-MIT-blue.svg)](https://github.com/Tharos-ux/wisp/blob/master/LICENCE)

<img align="right" src="https://github.com/Tharos-ux/wisp/blob/master/preview/WISP.png" alt="wisp logo" width="300"/>

This Python program is meant to determine to which taxa a bacteria is belonginig to, from long reads (>10.000bp), solely based upon alignment-free methods. As of now, the main focus is upon kmers proportions. It aims to do binning over a collection of samples, giving a probable class to each read.
Currently, five levels of taxa are implemented : **domain**, **phylum**, **group**, **order** and **family**.
Once a model finished at a given taxa level, it aims to do another iteration from previous results, excluding non-matching reference genomes.
This work was extensively tested both in baseline and leave-one-out scenarios. The two following figures are averaged over 103200 distinct reads, and depicts baseline accuracy with no fine tuning of parameters at family level. You may find more details in the linked publication below.


<p align="center">
  <img src="https://github.com/Tharos-ux/wisp/blob/master/preview/leaveoneout_supported_group.png" width="32%" />
  <img src="https://github.com/Tharos-ux/wisp/blob/master/preview/leaveoneout_supported_order.png" width="32%" /> 
  <img src="https://github.com/Tharos-ux/wisp/blob/master/preview/leaveoneout_supported_family.png" width="32%" />
</p>


The core functionnalities relies on a class probabiliy attribution to discriminate reads that might not be good indicators for our specie to be determined. As many other options, you can choose the ratio and the selection function to suit best your biological context.
A result for a sample read issued from a MinION lecture of Streptococcus may be found here : [MinION read](preview/sample_report.pdf)

**WISP is research software**. If you want to use it, please source the code. 

<p align="center">
  <img src="https://github.com/Tharos-ux/wisp/blob/master/preview/compdiff_3d_Actinobacteria.png" width="32%" />
  <img src="https://github.com/Tharos-ux/wisp/blob/master/preview/compdiff_3d_Fusobacteria.png" width="32%" />
  <img src="https://github.com/Tharos-ux/wisp/blob/master/preview/compdiff_3d_Proteobacteria.png" width="32%" /> 
</p>

It elaborates upon the notion of quantitative kmer signatures, and seeks to define patterns inside smaller unities (10.000 bp fragments). It was tested on errorless and noisy data, you may see a [study](preview/Bacterial_families_identification_from_long_reads___machine_learning_with_XGBoost_.pdf) where we go more in detail about it (originally written for [MLMG 2022](https://mlmg2022.github.io/)).
This tool is requiring some reference genomes, which it will index, to create a XGBoost model. Genomic fasta (.fna) files are the prefferd input style as of now. One can download custom genome dataset with NCBI accession numbers to create its own specific dataset and increase even more classifier accuracy. You may find here [accession files](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) that WISP can directly take as inputs to retrieve both genomes and taxonomy.

The remaining of this Readme will focus on WISP installation and setup, see the [wiki](https://github.com/Tharos-ux/wisp/wiki) for information 

# Installation

The easiest way to install and use WISP is through conda. A script to create the conda venv is included in this repo ; sipmly run the command and you will get a clean environnment for WISP, with all its dependencies :

```bash
bash env.sh
```

If you don't want to use conda, you can do a raw install of WISP dependencies : it requires a Python 3.10 version, and all packages listed in the *requirements.txt* file. Once your python envrionnement or installation is up, you can install quickly all dependencies with :

```bash
pip install -r requirements.txt
```

Also requires a Graphviz installation in order to perform trees rendering (currently disabled because of package issue with python 3.10.4).

# Submitting a job through conda

To easier the process of using the environment, disconecting from it and such, a bash file allows you to simply 

```bash
sbatch --cpus-per-task=8 --mem=50G wisp.sh ENV_PATH INTERFACE
```

+ ENVPATH is the path to your previously created conda environment
+ INTERFACE is the main call to a WISP subprogram, with its parameters

The next lines will describe interfaces you can call, with their parameters.

# Interface : WISP with wisp.py

It will iterate through your unknown genomes and create separate reports for each one, and build databases and models if they don't exist.
You can enforce database creation through this interface.

```bash
wisp.py [-b] [-t number]
```

**Optional arguments :**

+ -t, --multithreading > you can imput an int to define the number of parallel processes you want to spawn (recommanded : your number of cores)
+ -b, --build > tells WISP to build a full database instead of doing a prediction

# Interface : generating parameters file with parameters_init.py

You are free to change those parameters to better fit your current classification requirements.
Code and comments are self explanatory enough, so see this snippet of code for documentation.

```python
    sparkle: str = "small"
    # dict to be converted in .json file to create our parameters set
    params_job: dict = {
        # subreads size and max. of reads in sample
        'window_size': 10000,
        'sampling_objective': 500,
        # params for your database here [kmer_size, subsampling_depth, pattern]
        'domain_ref': [5, 50, [1, 1, 1, 1, 1]],
        'phylum_ref': [5, 100, [1, 1, 1, 1, 1]],
        'group_ref': [4, 100, [1, 1, 1, 1]],
        'order_ref': [4, 100, [1, 1, 1, 1]],
        'family_ref': [4, 100, [1, 1, 1, 1]],
        'merged_ref': [4, 50, [1, 1, 1, 1]],
        # params for your sample here [kmer_size, pattern]
        'domain_sample': [5, [1, 1, 1, 1, 1]],
        'phylum_sample': [5, [1, 1, 1, 1, 1]],
        'group_sample': [4, [1, 1, 1, 1]],
        'order_sample': [4, [1, 1, 1, 1]],
        'family_sample': [4, [1, 1, 1, 1]],
        'merged_sample': [4, [1, 1, 1, 1]],
        # 'input' : location of reference genomes
        'input_train': f"genomes/train_{sparkle}/",
        # 'input_unk' : location of unk genomes
        'input_unk': f"genomes/unk_{sparkle}/",
        # 'output' : output for database
        'database_output': "data/",
        'reports_output': f"output/{sparkle}/",
        # parameters for exploration and algorithm
        'threshold': 0.10,
        'nb_boosts': 10,
        'tree_depth': 10,
        # parameters regarding results
        'test_mode': 'min_set',  # 'no_test', 'min_set', 'verbose'
        # parameter for read selection, significance for softprob
        'reads_th': 0.1,
        'selection_mode': 'delta_mean',  # 'min_max','delta_mean','delta_sum'
        # force rebuilding full model, when you re-use a database but you changed model parameters
        'force_model_rebuild': True,  # never set true in multithread mode
        # tells the software if should consider both orientations or assume it is 3' -> 5' and computes canonical kmers
        'single_way': True,
        'targeted_level': 'family',  # domain, phylum, group, order, family
        'levels_list': ['domain', 'phylum', 'group', 'order', 'family'],
        'abundance_threshold': 0.25,
        # to fetch WISP genomes from refseq
        'email': 'XXXXXXXXXXXXXXX@XXXXXXX.XX',
        'annotate_path': 'genomes/to_annotate',
        'accession_numbers': 'genomes/assembly_summary.txt',
        # name for database
        'db_name': sparkle,
        'prefix_job': sparkle,
        'log_file': f'LOG_wisp_{sparkle}'
    }
```

You can modify the raw .json file, or update the dict contained in script *wisp_lib/parameters_init.py* and re-generate the file with :

```bash
wisp_lib/parameters_init.py filename
```

**Optional arguments :**

+ -t, --multithreading > you can imput an int to define the number of parallel processes you want to spawn (recommanded : your number of cores)

# Results

Results are currently exported into a set of plots, a tex file and a json file, in order to easily overview data or fetch it to another program.

**How to interpret a report**

For each classification level, you get a set of infos. You may get even more with `'full_test_set': verbose` but compute time will be higher. The first row tells the global classification result, based upon the reads selection method you choose.
At the end of the report, you may find the full parameters set you used on this particular sample, for archiving purposes.
The tree displays information on which ways were investigated, matching the threshold you defined. Nodes are classes and edges are percentages of attributed reads on those classes from previous node.

Then, you get the information by level : in base mode, you have one confusion matrix, displaying the results from a test set against the trained model ; numbers are row percentages, showing for each actual class the percentage of correctly (or not) predicted reads. Note that this info is calculated upon your selection parameters, so you can investigate which parameters works better for your database, but you should look out for overfitting!

The next graph tells you about reads attribution : once calculations are done, reads are attributed to a certain class.
The dashed line displays the investigation limit you set (reads percentage to consider a class)

Next up, if you choose advanced plotting options, you get an additionnal set of graphs.
Firstly, you have 'Mean and standard deviation' which tells you more about the fidelity of successive boostings.

Secondly, you get a plotting of the 15 most important features for the split, to get an extra bit more info upon the signature concept which this algorithm relies on.

Lastly, you get additional plots that are not included inside the .tex file but are freely consultable from the output folder :

+ Reads probabilities repartition
+ Reads selection across functions