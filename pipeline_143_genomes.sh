#!/bin/sh
#SBATCH --job-name=pipeline_143_genomes
#SBATCH --output=LOG_job.log
. /local/env/envconda.sh

# create mandatory dirs if not exist
mkdir -p genomes/unk
mkdir -p genomes/train
mkdir -p data
mkdir -p output

# use conda env
conda activate "/home/genouest/genscale/sdubois/wisp-env"

for FILE in genomes/train_small/*
do
    mv "genomes/train_small/"$(basename $FILE) "genomes/unk_small/"
    rm -r data/tmp_database
    python wisp.py -b -t 8
    python wisp.py
    mv "genomes/unk_small/"$(basename $FILE) "genomes/train_small/"
done

conda deactivate