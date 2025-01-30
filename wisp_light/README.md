# Build et Predict
depuis genouest ou sur laptop
## env
```
srun --time 00-10:00:00 --pty bash   #depuis genouest
. ~/envtaxo2/bin/activate
cd code/wisp/workspace
```

## Build 
sur toute la base refseq genouest , un peu long, à lancer depuis submit_main_build.sh avec suffisament de RAM

> python main.py build refseq /groups/microtaxo/data/refseq_with_taxo/

# Predict
> python main.py predict refseq /groups/microtaxo/data/refseq_with_taxo/ /home/genouest/cnrs_umr6074/hcourtei/out_refseq

from laptop
> python main.py predict refseq /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo /home/hcourtei/Projects/MicroTaxo/codes/data/out_refseq


TODO

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
