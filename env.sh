#!/bin/sh
#SBATCH --job-name=wisp
. /local/env/envconda.sh

# init env
conda create -p pwd"/wisp_env" python=3.10

conda activate pwd"/wisp_env"
# installing packages
python -m pip install -r requirements.txt

conda deactivate