#!/bin/sh
#SBATCH --job-name=wisp
#SBATCH --output=logs/wisp-%j.log
#SBATCH --mail-user=siegfried.dubois@inria.fr
#SBATCH --mail-type=END,FAIL

. /local/env/envconda.sh

# create mandatory dirs if not exist
mkdir -p genomes/unk
mkdir -p genomes/train
mkdir -p data
mkdir -p output
mkdir -p logs

# use conda env
WD=$(pwd)
conda activate $WD"/wisp_env"

python $1

conda deactivate

# example use:
# sbatch --cpus-per-task=8 --mem=80G wisp.sh "wisp.py parameters_files/megadb -b"
# sbatch wisp.sh "utilities.py extract_genomes 3 /scratch/sdubois/refseq /scratch/sdubois/lb_minion/train_mega"