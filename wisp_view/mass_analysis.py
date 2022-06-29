"Aim of this module is to do comparaisons between sets of parameters, assuming the same series were run multiple times with different parameters"

from collections import Counter
from os import listdir, system
from json import load
from numpy import outer
from pandas import DataFrame
from random import shuffle
from argparse import ArgumentParser

FOLDER_LIST: list[str] = [
    '10p_010dmean_143genomes',
    '10p_010dmean_143genomes_destr50percent',
    '10p_010dmean_400genomes',
    '10p_010dmean_400genomes_destr5percent',
    '10p_010dmean_400genomes_destr50percent'
]


def compare():
    df = DataFrame()
    for subfolder in FOLDER_LIST:
        temp_list = []
        for report_folder in listdir(f"output/{subfolder}/"):
            try:
                json_report = [file for file in listdir(
                    f"output/{subfolder}/{report_folder}/") if '.json' in file][0]
                with open(f"output/{subfolder}/{report_folder}/{json_report}") as json_reader:
                    temp_dict = load(json_reader)
                    temp_list.append(
                        f"{temp_dict['phylum']}, {temp_dict['group']}, {temp_dict['order']}, {temp_dict['family']}")
            except:
                temp_list.append('-')  # lacking report file
        df[subfolder] = temp_list
        df.index = listdir(f"output/{subfolder}/")
    df.to_csv("output/mass_report.csv")


PATH_TRAIN: str = "/udd/sidubois/Stage/output/onevsall_sampled_v1/"
PATH_SAMPLED: str = "/udd/sidubois/Stage/output/onevsall_sampled_v1/"
LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']


def number_of_classes(t: int, input_path: str, output_path: str):
    """Computes the number of classes at each level
    """
    t = int(t)
    extract_families_with_enough_representatives(t, input_path, output_path)
    for direct in [input_path, output_path]:
        all_genomes: list[list] = [
            [e.split('_')[i] for e in listdir(direct)] for i in range(len(LEVELS))]
        for i, level in enumerate(all_genomes):
            counts = Counter(level)
            print(
                f"{LEVELS[i]} : {len(counts)} classes, with {sum([1 for _,v in counts.items() if v>=t])} equal or above {t} representatives")


def extract_families_with_enough_representatives(t: int, input_path: str, output_path: str):
    t = int(t)
    all_genomes: list = [e.split('_')[4] for e in listdir(input_path)]
    print(f"Database contains {len(all_genomes)} reference genomes")
    selected = [f for f, v in Counter(all_genomes).items() if v >= t]
    print(f"Len: {len(selected)} -> {selected}")
    if t != 0:
        extraction(selected, input_path, output_path)


def extraction(selected: list, input_path: str, output_path: str):
    all_genomes_full: list = [e for e in listdir(PATH_TRAIN)]
    shuffle(all_genomes_full)
    for enf in selected:
        list_selected = [g for g in all_genomes_full if enf == g.split('_')[
            4]][0:3]
        for file in list_selected:
            system(f"cp {input_path}/{file} {output_path}/{file}")
