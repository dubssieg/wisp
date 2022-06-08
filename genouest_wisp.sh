#!/bin/sh
#SBATCH --job-name=wisp
. /local/env/envconda.sh

conda activate $1

##########################
#                        #
#    ARGS YOU CAN USE    #
#                        #
##########################

# calling build
# python wisp.py -b -t 8

# calling dl
# python utilities.py

# calling prediction
# python wisp.py -t 8

$2

conda deactivate

# example use:
# sbatch --cpus-per-task=8 --mem=50G genouest_wisp.sh "/home/genouest/genscale/sdubois/wisp-env" "wisp.py -b -t 8"