#!/bin/sh
#SBATCH --job-name=wisp
. /local/env/envconda.sh

# init env
conda create wisp_env python

conda activate wisp_env

# installing packages
pip install -r requirements.txt

conda deactivate