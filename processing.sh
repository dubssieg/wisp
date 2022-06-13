#!/bin/bash

for FILE in genomes/train_small/*
do
    mv "genomes/train_small/"$(basename $FILE) "genomes/unk_small/"
    rm -r data/tmp_database
    python wisp.py -b -t 8
    python wisp.py
    mv "genomes/unk_small/"$(basename $FILE) "genomes/train_small/"
done