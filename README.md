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

# Quick setup : see [Wiki](https://github.com/Tharos-ux/wisp/wiki/User-guide) for detailed instructions

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