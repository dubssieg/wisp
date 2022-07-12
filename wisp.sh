#!/bin/sh
#SBATCH --job-name=wisp
#SBATCH --output=logs/wisp-%j.log
#SBATCH --mail-user=siegfried.dubois@irisa.fr
#SBATCH --mail-type=END,FAIL

. /local/env/envconda.sh

# create mandatory dirs if not exist
mkdir -p genomes/unk
mkdir -p genomes/train
mkdir -p data
mkdir -p output

# use conda env
conda activate $1

##########################
#                        #
#    ARGS YOU CAN USE    #
#                        #
##########################

# calling build
# wisp.py -b -t 8

# calling dl
# utilities.py

# calling prediction
# wisp.py -t 8

python $2

conda deactivate

# example use:
# sbatch --cpus-per-task=8 --mem=80G wisp.sh "/home/genouest/genscale/sdubois/wisp-env" "wisp.py parameters_files/megadb -b"