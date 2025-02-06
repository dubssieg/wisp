# train xgboost on refseq

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



# Old wisp command
## Build
sur toute la base refseq genouest , un peu long, à lancer depuis submit_main_build.sh avec suffisament de RAM

> python main.py build refseq /groups/microtaxo/data/refseq_with_taxo/

## Predict
> python main.py predict refseq /groups/microtaxo/data/refseq_with_taxo/ /home/genouest/cnrs_umr6074/hcourtei/out_refseq

from laptop
> python main.py predict refseq /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo /home/hcourtei/Projects/MicroTaxo/codes/data/out_refseq


# TODO

## examen base
- version refseq  :  Release 227 November 4, 2024.
- filtrage refseq stat pattern 6 level domain -> specie
- faire stat sur refseq, prendre un representant par espèces , nb genome par famille 
- retenir famille si au moins 10 représentant
- nb famille dans ce cas

## split train/val
seed pour random
entrée : 1 liste de path vers .fna
Sur 10 représentant dans une famille , en mettre 1 dans en val
sortie 1 liste en train / 1 liste en val
