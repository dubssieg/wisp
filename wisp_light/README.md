# train xgboost on refseq
auteur package wisp_light: Hermann Courteille (PNRIA)

projet: micro taxo:
lancement : 25 novembre 2024
mi projet : 25 février
fin prévu : 25 mai

- version refseq  :  Release 227 November 4, 2024.


## make env
depuis le répertoire wisp_light

```
srun --time 00-10:00:00 --pty bash 
. /local/env/envpython-3.11.9.sh
python3.11 -m venv ~/envtaxo2
source ~/envtaxo2/bin/activate
pip install -r requirements.txt 
```

# Train and evaluate xgboost models


data are by default:
`datadir=/projects/microtaxo/data/refseq_with_taxo_merged`
See and edit : `params.yaml`

Train
```
srun --time 00-10:00:00 --mem=20G --cpus-per-task=8 --pty bash #depuis genouest
source ~/envtaxo2/bin/activate
cd ~/codes/wisp/wisp_light/training
python train_val.py 
```

Restart training from existing json database:

> python train_val.py --db_json  /home/genouest/cnrs_umr6074/hcourtei/codes/wisp/exp/model_base_02_05_16_55/databases.json

Logs and result are in exp_rootdir by default
`~/codes/wisp/exp`

1. Interactive session with srun above :
To prevent ssh break, you can use tmux on genouest see https://help.genouest.org/usage/slurm/#long-running-interactive-jobs-srun

2. Sbatch , fix parameter in .sh , params.yaml or train_val.py, then 
`sbatch submit_main_build.sh
`
# See results 

from compute  <node_name>  in genouest:

>mlflow ui--port 8123 --backend-store-uri /projects/microtaxo/exp_refseq/mlruns

from local laptop

>ssh -A -t -t hcourtei@genossh.genouest.org -L 8123:localhost:8123 ssh <node_name> -L 8123:localhost:8123



# conda env with glibc >1.28 

> 2.1 it/s

`
conda install -y gcc_linux-64 gxx_linux-64 -c conda-forge
pip install xgboost --no-binary :all:
`
```
. /local/env/envconda.sh
conda activate py311_env
source ~/.bashrc
```

```.bashrc

export PATH=$CONDA_PREFIX/libexec/gcc/x86_64-conda-linux-gnu/14.2.0:$PATH
export CC=$CONDA_PREFIX/libexec/gcc/x86_64-conda-linux-gnu/14.2.0/gcc
export CXX=$CONDA_PREFIX/libexec/gcc/x86_64-conda-linux-gnu/14.2.0/g++
```

conda install -c nvidia cudatoolkit=11.8.0

for nvcc

export PATH=/usr/local/cuda-12.3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-12.3/lib64:$LD_LIBRARY_PATH
nvcc --version 

## test avec gpu

srun --time 00-01:00:00 --mem=20G --gpus 1 -p gpu --pty bash



# Old wisp command
## Build
sur toute la base refseq genouest , un peu long, à lancer depuis submit_main_build.sh avec suffisament de RAM

> python main.py build refseq /groups/microtaxo/data/refseq_with_taxo/

## Predict
> python main.py predict refseq /groups/microtaxo/data/refseq_with_taxo/ /home/genouest/cnrs_umr6074/hcourtei/out_refseq

from laptop
> python main.py predict refseq /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo /home/hcourtei/Projects/MicroTaxo/codes/data/out_refseq


# TODO
