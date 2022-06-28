#!/bin/sh
#SBATCH --job-name=wisp_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --output=LOG_pipeline.log
#SBATCH --mail-user=siegfried.dubois@irisa.fr
#SBATCH --mail-type=END,FAIL

# sourcing conda
. /local/env/envconda.sh

# variables
ENV="/home/genouest/genscale/sdubois/wisp-env"
PARAMETERS="parameters_files/"$1".json"
INPUT_MINION="/scratch/sdubois/Lactobacillales_MinION/pass"
OUTPUT_MINION="genomes/unk_"$1
REF_GENOMES="genomes/train_"$1

# create mandatory dirs if not exist
mkdir -p $OUTPUT_MINION
mkdir -p $REF_GENOMES
mkdir -p "data"
mkdir -p "output/"$1

# use conda env
conda activate $ENV

    # pre-processing .fastq MinION reads
    python utilities.py 'clean_minion' $INPUT_MINION $OUTPUT_MINION

    # downloading set of reference genomes
    python utilities.py 'summary_to_dl' $SUMMARY_FILE $REF_GENOMES

    # cleaning somewhat dirty NCBI classification
    python utilities.py 'clean_rename' $REF_GENOMES

    # calling for db creation
    python wisp.py $PARAMETERS -b

    # calling for analysis
    python wisp.py $PARAMETERS -t 8

    # plotting facilities for visualisation
    # to be added

# exiting conda env
conda deactivate