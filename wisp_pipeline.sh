#!/bin/sh
#SBATCH --job-name=wisp_pipeline
#SBATCH --output=LOG_pipeline.log
#SBATCH --mail-user=siegfried.dubois@irisa.fr
#SBATCH --mail-type=END,FAIL
. /local/env/envconda.sh

# variables

INPUT_MINION="/scratch/sdubois/Lactobacillales_MinION/pass"

# $1 must be something like "subsampled","full","143genomes"

mkdir -p "genomes/MinION_reads"

# create mandatory dirs if not exist
#mkdir -p "genomes/unk_"$1
#mkdir -p "genomes/train_"$1
#mkdir -p "data"
#mkdir -p "output/"$1


# use conda env
conda activate "/home/genouest/genscale/sdubois/wisp-env"

python utilities.py 'clean_minion' $INPUT_MINION

# python wisp.py "wisp_params_"$1".json" -b -t 8

#for FILE in "genomes/train_"$1"/"*
#do
#    mv "genomes/train_"$1"/"$(basename $FILE) "genomes/unk_"$1"/"
#    python wisp.py "wisp_params_"$1".json"
#    mv "genomes/unk_"$1"/"$(basename $FILE) "genomes/train_"$1"/"
#done

conda deactivate