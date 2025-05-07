

**projet**: microtaxo

**accompagnement PNRIA**: du 25 novembre 2024 au 25 mai 2025


**package wisp_light**: Hermann Courteille (PNRIA)

# I. Set Environment
## conda
```
conda env create -f micro_env.yml
conda activate micro_env
```
## virtual env
sur genouest, obligatoirement sur un noeud calcul
```
srun --time 00-10:00:00 --pty bash 
. /local/env/envpython-3.11.9.sh

python3.11 -m venv ~/envtaxo
source ~/envtaxo/bin/activate
pip install -r requirements.txt 
```



# II. Build refseq dataset
## a. Télécharger et dézipper tous les fichiers listés dans le .tsv 
- obtenir le summary.txt 
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -o {output_dir}
```
ou via l'interface web 
https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2&reference_only=true
- wisp_light/build_dataset/refseq/reference_genome_summary.tsv
**version refseq**:  Release 227 November 4, 2024.

depuis build_refseq/refseq
```
python import_refseq.py --csv_file reference_genome_summary.tsv --output_dir /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data --num_workers 8
```
### b. Obtenir les taxonomies à partir des taxid de chaque génome 

```
python get_all_taxo_from_NCBI.py --input reference_genome_summary.tsv --taxid_column taxid --batch_size 10
```
2 tableaux tsv sont générés : 
- le 1er avec toutes les taxonomies présentes  ['phylum', 'class', 'order', 'family'] dans "reference_genome_summary_complete_taxo.tsv"
- le 2eme avec des taxonomies incompletes dans "reference_genome_summary_incomplete_taxo.tsv"

### c. association dans l'itérateur
```python
from wisp_light.dataset.refSeqDataset import RefSeqDataset
datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data"
dataset = RefSeqDataset("reference_genome_summary_complete_taxo.tsv", datadir)
```
```
14:18 - refSeqDataset.py - INFO - nb files in datadir: 451 restrict to 448 with labels in index 
14:18 - refSeqDataset.py - INFO - Splitted dataset nb 448 into train :403 val: 45
 0, genome /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data/GCF_002865995.1_ASM286599v1_genomic.fna
taxo_dict {'phylum': 'Bacillota', 'class': 'Clostridia', 'order': 'Peptostreptococcales', 'family': 'Peptostreptococcaceae'}
 1, genome /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data/GCF_016028775.1_ASM1602877v1_genomic.fna
taxo_dict {'phylum': 'Pseudomonadota', 'class': 'Gammaproteobacteria', 'order': 'Moraxellales', 'family': 'Moraxellaceae'}
```
# III. Train and evaluate xgboost models
On part de tous les génomes de références de refseq.
Cette base est découpée en train/val avec la fonction sklearn.model_selection.train_test_split, la seed donné
dans params.yaml afin de pouvoir être reproduit.

L'entrainement et la validation se fait en 3 temps:
1. la construction d'une base de données regroupant tous les comptages de kmer sur le train
2. l'entrainement des modèles de façon hiérachique
3. la validation sur le jeu de validation

Tous les paramètres de pre-processing, de xgboost ... sont dans : `wisp_light/training/params.yaml`

Commande à partir de wips_light
en local : 
```
export PYTHONPATH=/home/hcourtei/Projects/MicroTaxo/codes

 python wisp_light/training/train_val.py  \
  --exp_name test_laptop  \
  --datadir /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data \
  --params_file training/params.yaml
```
sur genouest, les données sont sur /projects/microtaxo/data
.
├── AllTheBacteria
├── refseq_complete_genome
└── refseq_reference_genome


```
srun --time 00-10:00:00 --mem=20G --cpus-per-task=8 --pty bash #depuis genouest
source ~/envtaxo2/bin/activate
cd ~/codes/wisp

python wisp_light/training/train_val.py \
  --exp_name test_laptop  \
  --datadir /projects/microtaxo/data/refseq_reference_genome \
  --params_file training/params.yaml
  ```

Restart training from existing json database:

> python train_val.py --db_json  /home/genouest/cnrs_umr6074/hcourtei/codes/wisp/exp/model_base_02_05_16_55/databases.json

Logs and result are in exp_rootdir by default
`~/codes/wisp/exp`

1. Interactive session with srun above :
To prevent ssh break, you can use tmux on genouest see https://help.genouest.org/usage/slurm/#long-running-interactive-jobs-srun

2. Sbatch , fix parameter in .sh , params.yaml or train_val.py, then 
`sbatch submit_main_build.sh`

# IV. See results 

depuis un <noeud> de calcul  in genouest:
```
tmux # pour avoir une session détachée
srun --pty --time=08:00:00 bash
. ~/envtaxo2/bin/activate

mlflow ui --port 8123 --backend-store-uri /projects/microtaxo/exp_refseq/mlruns
```
from local laptop
```
ssh -A -t -t hcourtei@genossh.genouest.org -L 8123:localhost:8123 ssh <noeud> -L 8123:localhost:8123
ls /projects/microtaxo/exp_refseq/
```

print(list(counters[0].items())[:10])


partition avec disque plus rapide

srun --cpus-per-task=20 -p genscale -w cl1n027 --mem 40600 --pty bash

- cl1n026 (24 Xeon(R) CPU E5-2640 0 @ 2.50GHz)
- cl1n027 (40  Xeon(R) CPU E5-2660 v3 @ 2.60GHz)
- cl1n028 (40  Xeon(R) CPU E5-2660 v3 @ 2.60GHz)

# DIVERS


## conda env with glibc >1.28 

> 2.1 it/s

```
conda install -y gcc_linux-64 gxx_linux-64 -c conda-forge
pip install xgboost --no-binary :all:
```

```
. /local/env/envconda.sh
conda activate py311_env
source ~/.bashrc
```

```
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
