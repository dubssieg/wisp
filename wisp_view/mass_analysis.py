"Aim of this module is to do comparaisons between sets of parameters, assuming the same series were run multiple times with different parameters"

from collections import Counter
from os import listdir, system
from json import load
from pandas import DataFrame
from random import shuffle

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
PATH_SAMPLED: str = "genomes/sampled/"
LEVELS: list[str] = ['domain', 'phylum', 'group', 'order', 'family']


def number_of_classes(t: int):
    """Computes the number of classes at each level
    """
    extract_families_with_enough_representatives(t)
    for dir in [PATH_TRAIN, PATH_SAMPLED]:
        all_genomes: list[list] = [
            [e.split('_')[i] for e in listdir(dir)] for i in range(len(LEVELS))]
        for i, level in enumerate(all_genomes):
            counts = Counter(level)
            print(
                f"{LEVELS[i]} : {len(counts)} classes, with {sum([1 for _,v in counts.items() if v>=t])} equal or above {t} representatives")


def extract_families_with_enough_representatives(t: int):
    all_genomes: list = [e.split('_')[4] for e in listdir(PATH_TRAIN)]
    print(f"Database contains {len(all_genomes)} reference genomes")
    selected = [f for f, v in Counter(all_genomes).items() if v >= t]
    print(f"Len: {len(selected)} -> {selected}")
    # extraction(selected)


def extraction(selected):
    all_genomes_full: list = [e for e in listdir(PATH_TRAIN)]
    shuffle(all_genomes_full)
    for f in selected:
        list_selected = [g for g in all_genomes_full if f == g.split('_')[
            4]][0:3]
        for e in list_selected:
            system(f"cp genomes/train/{e} genomes/sampled/{e}")


number_of_classes(3)
