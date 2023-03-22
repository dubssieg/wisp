"Comparaisons between sets of parameters, assuming the same series were run multiple times"

from collections import Counter
from os import listdir, system
from json import load
from random import shuffle
from pandas import DataFrame


def compare(folder_list: list) -> None:
    """Aggregates comparaison outputs of folders inside folder

    Args:
        folder_list (list): all folders we want to compare
    """
    df = DataFrame()
    for subfolder in folder_list:
        temp_list = []
        for report_folder in listdir(f"output/{subfolder}/"):
            try:
                json_report = [file for file in listdir(
                    f"output/{subfolder}/{report_folder}/") if '.json' in file][0]
                with open(f"output/{subfolder}/{report_folder}/{json_report}", encoding="utf-8") as json_reader:
                    temp_dict = load(json_reader)
                    temp_list.append(
                        f"{temp_dict['phylum']}, {temp_dict['group']}, {temp_dict['order']}, {temp_dict['family']}")
            except Exception:
                temp_list.append('-')  # lacking report file
        df[subfolder] = temp_list
        df.index = listdir(f"output/{subfolder}/")
    df.to_csv("output/mass_report.csv")


def number_of_classes(targeted: int, input_path: str, output_path: str) -> None:
    """Computes the number of classes at each level

    Args:
        targeted (int): number of classes we want to extract
        input_path (str): path where genomes are
        output_path (str): path where extracted genomes will go
    """
    levels: list[str] = ['domain', 'phylum', 'group', 'order', 'family']
    targeted = int(targeted)
    extract_families_with_enough_representatives(
        targeted, input_path, output_path)
    for direct in [input_path, output_path]:
        all_genomes: list[list] = [
            [e.split('_')[i] for e in listdir(direct)] for i in range(len(levels))]
        for i, level in enumerate(all_genomes):
            counts = Counter(level)
            print(
                f"{levels[i]} : {len(counts)} classes, with {sum([1 for _,v in counts.items() if v>=targeted])} equal or above {targeted} representatives")


def extract_families_with_enough_representatives(targeted: int, input_path: str, output_path: str) -> None:
    """Extract families matching threshold

    Args:
        targeted (int): number of classes we want to extract
        input_path (str): path where genomes are
        output_path (str): path where extracted genomes will go
    """
    targeted = int(targeted)
    all_genomes: list = [e.split('_')[4] for e in listdir(input_path)]
    print(f"Database contains {len(all_genomes)} reference genomes")
    selected = [f for f, v in Counter(all_genomes).items() if v >= targeted]
    print(f"Len: {len(selected)} -> {selected}")
    if targeted != 0:
        extraction(selected, input_path, output_path)


def extraction(selected: list, input_path: str, output_path: str) -> None:
    """Extract genomes that are listed into another folder

    Args:
        selected (list): selected genomes
        input_path (str): path where to take them
        output_path (str): path where to copy them
    """
    all_genomes_full: list = [e for e in listdir(input_path)]
    shuffle(all_genomes_full)
    for enf in selected:
        list_selected = [g for g in all_genomes_full if enf == g.split('_')[
            4]][0:3]
        for file in list_selected:
            system(f"cp {input_path}/{file} {output_path}/{file}")
