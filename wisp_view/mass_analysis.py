# Aim of this module is to do comparaisons between sets of parameters, assuming the same series were run multiple times with different parameters
from os import listdir
from json import load
import pandas as pd

list_of_folders = ['purge_reads_0,00_8_boosts',
                   'purge_reads_0,00_10_boosts',
                   'purge_reads_0,50_8_boosts',
                   'purge_reads_0,50_10_boosts',
                   'purge_reads_0,25_8_boosts',
                   'purge_reads_0,25_10_boosts',
                   'purge_reads_0,25_12_boosts']


def compare():
    df = pd.DataFrame()
    for subfolder in list_of_folders:
        temp_list = []
        for report_folder in listdir(f"output/{subfolder}/"):
            json_report = [file for file in listdir(
                f"output/{subfolder}/{report_folder}/") if '.json' in file][0]
            with open(f"output/{subfolder}/{report_folder}/{json_report}") as json_reader:
                temp_dict = load(json_reader)
                temp_list.append(
                    f"{temp_dict['phylum']}, {temp_dict['group']}, {temp_dict['order']}, {temp_dict['family']}")
        df[subfolder] = temp_list
        df.index = listdir(f"output/{subfolder}/")
    df.to_csv("output/mass_report.csv")


compare()
