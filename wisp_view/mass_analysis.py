"Aim of this module is to do comparaisons between sets of parameters, assuming the same series were run multiple times with different parameters"

from collections import Counter
from os import listdir
from json import load
from pandas import DataFrame

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


PATH_TRAIN: str = "genomes/train/"
LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']


def number_of_classes():
    all_genomes: list[list] = [
        [e.split('_')[i] for e in listdir(PATH_TRAIN)] for i in range(len(LEVELS))]
    for i, level in enumerate(all_genomes):
        print(f"{LEVELS[i]} : {len(Counter(level))} classes")


# compare()

number_of_classes()
