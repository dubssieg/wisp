#!/bin/sh
#SBATCH --job-name=wisp
. /local/env/envconda.sh

# init env
conda create -p wisp_env python=3.10

conda activate wisp_env

# installing packages
pip install -r requirements.txt

conda deactivate