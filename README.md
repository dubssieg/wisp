<p align="center">
[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![](https://img.shields.io/badge/documentation-unfinished-orange.svg)]()
[![](https://img.shields.io/badge/comments-finished-green.svg)]()
[![](https://img.shields.io/badge/build-stable-green.svg)]()
</p>

# WISP : binning of bacterial long-read signatures

<img align="right" src=https://github.com/Tharos-ux/wisp/blob/master/preview/WISP.png alt="wisp logo" width="300"/>

This Python program is meant to determine to which taxa a bacteria is belonginig to, from long reads (>10.000bp), solely based upon alignment-free methods. As of now, the main focus is upon kmers proportions. It aims to do binning over a collection of samples, giving possible class to each read.
It is a rebranching of a [master 1 internship project](https://github.com/Tharos-ux/wisp/tree/master), done on free time. Concept has not evolved since, but code was redesign for better comprehension.

Currently, five levels of taxa are implemented : **domain**, **phylum**, **group**, **order** and **family**.
Once a model finished at a given taxa level, it aims to do another iteration from previous results, excluding non-matching reference genomes.

The core functionnalities relies on a class probabiliy attribution to discriminate reads that might not be good indicators for our specie to be determined. As many other options, you can choose the ratio and the selection function to suit best your biological context.

**WISP is research software**. If you want to use it, please source the code. 

## Installing software

```bash
git clone -b v0.1.0 --single-branch git@github.com:Tharos-ux/wisp.git
cd wisp/
python -m pip install . --quiet
```

## Usage

```bash
usage: wisp [-h] [-l] {build,predict} ...

Bacteria family identification tool.

Subcommands:
  {build,predict}  Available subcommands
    build          Creates the database from the specified set of files.
    predict        Creates the samples and evaluates them.

Global Arguments:
  -h, --help       show this help message and exit
  -l, --locals     Display locals on error.
```

The command `build` allows to create models from a set of reference genomes.

```bash
usage: wisp build [-h] [-p PARAMETERS] database_name input_folder

positional arguments:
  database_name         Name for database
  input_folder          Input folder containig reference genomes

options:
  -h, --help            show this help message and exit
  -p PARAMETERS, --parameters PARAMETERS
                        Specifies a parameter file
```

The command `predict` offers to predict taxonomy of sample from computed models.

```bash
usage: wisp predict [-h] [-p PARAMETERS] database_name input_folder output_folder

positional arguments:
  database_name         Name for database
  input_folder          Input folder containig unknown genomes
  output_folder         Input folder containig reference genomes

options:
  -h, --help            show this help message and exit
  -p PARAMETERS, --parameters PARAMETERS
                        Specifies a parameter file
```

## Project architecture

All code about binning is in the `workspace` folder.
- `main.py` is the main loop and argument parser
- `create_database.py` contains function to index the reference genomes
- `create_model.py` contains the functions to create XGboost models from the index
- `create_sample.py` contains the functions to create the dataset for the reads we want to predict
- `create_prediction.py` contains the functions to make prediction on the sample dataset with models

Two scripts come along, in the `scripts` folder.
- `download_refseq.py` downloads, from a refseq assembly file, the representative genomes, and annotates them by thier classification (NCBI taxonomy)
- `visualize_output.py` renders a html file with graphs from a .json, output of the `wisp predict` command

URL to [refseq assembly file for **bacteria**](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt)

## Core idea

Given a set of reference genomes, annotated with their taxonomy, this program samples a set of 10kb lectures in each reference genome. Each sample is fractured in k-mers, which are counted : those counts are the features for the model. We train the model on many samples, in order to try to extract a k-mer signature for each specie in our reference genomes.
Then, when we want to predict, we apply the same treatment to our lectures, and we try to compute a list of classes that can be associated to the unknown sample. Once the upper classification level is determined, we move on to a lower one, until we reach family. If the confidence score is high enough, we may explore multiple branches of the taxonomy.

## Tasks

- [x] Proof-of-concept
- [x] Minimal docs + docstings
- [x] Add multithreading on database build (#2fab1e3)
- [x] Add multithreading on prediction (#9dbf868)
- [x] Report system (#c07c7af)
- [x] Build files + new argument parser (#9a6a118)
- [ ] Put the remaining parameters in code in the parameters file
- [ ] Test coverage
- [ ] Validation of parameters + definition of default parameters in parameters file
- [ ] Explore if using the full DNA alphabet (implemented) is better than solely A, T, C and G
- [ ] Validate scalability (database build holds with 14k reference genomes/prediction with 1M+ reads)
- [ ] Find ways to reject less reads (assembly?) because all reads inferior to length threshold are rejected
- [ ] Otherwise, find appropriate parameters for shorter reads, and validate those
- [ ] Increase speed by computing predictions by batches instead of one-by-one
- [ ] Validate quality of results with mock communities, damaged mock communities, then real data