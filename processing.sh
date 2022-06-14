#!/bin/bash

for FILE in "genomes/train_"$1/*
do
    mv "genomes/train_"$1"/"$(basename $FILE) "genomes/unk_small/"
    echo "Processing "$(basename $FILE)
    python wisp.py "wisp_params_"$1".json" -e $(basename $FILE) -v
    mv "genomes/unk_"$1"/"$(basename $FILE) "genomes/train_small/"
done