#!/bin/sh
#SBATCH --job-name=wisp
. /local/env/envconda.sh

WD=$(pwd)

mkdir -p "data"
mkdir -p "output/"
mkdir -p "logs"
mkdir -p "genomes/train_default"
mkdir -p "genomes/unk_default"
mkdir -p "parameters_files"

# init env
conda create -p $WD"/wisp_env" python=3.10

conda activate $WD"/wisp_env"
# installing packages
python -m pip install -r requirements.txt

conda deactivate