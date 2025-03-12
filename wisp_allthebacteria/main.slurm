#!/bin/bash
#SBATCH --job-name=clx_create_db
#SBATCH --chdir=/home/genouest/cnrs_umr6074/cleroux2/microtaxo/wisp/wisp_allthebacteria
#SBATCH --output=slurm/%j.txt
#SBATCH --error=slurm/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=7-00:00:00 
##SBATCH --mem-per-cpu=20000
#SBATCH --mem=30G
#SBATCH --mail-user=cyrille.leroux@irisa.fr
#SBATCH --mail-type=END,FAIL


source ~/.bashrc
conda activate wisp

python main.py --create-db
