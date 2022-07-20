#!/bin/sh
#SBATCH --job-name=wisp_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --output=logs/pipeline-%j.log
#SBATCH --mail-user=siegfried.dubois@irisa.fr
#SBATCH --mail-type=END,FAIL

# sourcing conda
. /local/env/envconda.sh

# variables
WD = pwd
ENV=$WD"/wisp-env"
PARAMETERS="parameters_files/"$1".json"
INPUT_MINION="genomes/Lactobacillales_MinION/pass"
OUTPUT_MINION="genomes/unk_"$1
TEMP_GENOMES="genomes/temp_"$1
REF_GENOMES="genomes/train_"$1
NUMBER_RELATIVES=3

# create mandatory dirs if not exist
mkdir -p $OUTPUT_MINION
mkdir -p $REF_GENOMES
mkdir -p $TEMP_GENOMES
mkdir -p "data"
mkdir -p "output/"$1
mkdir -p "logs"

# use conda env
conda activate $ENV

    # pre-processing .fastq MinION reads
    python utilities.py 'clean_minion' $INPUT_MINION $OUTPUT_MINION

    # downloading set of reference genomes
    python utilities.py 'summary_to_dl' $SUMMARY_FILE $TEMP_GENOMES

    # selecting only ref genomes
    python utilities.py 'extract_genomes' $NUMBER_RELATIVES $TEMP_GENOMES $REF_GENOMES

    # cleaning somewhat dirty NCBI classification
    python utilities.py 'clean_rename' $REF_GENOMES

    # calling for db creation
    python wisp.py $PARAMETERS -b

    # calling for analysis
    python wisp.py $PARAMETERS

# exiting conda env
conda deactivate