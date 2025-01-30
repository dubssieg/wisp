import os
import pandas as pd

import natsort
import matplotlib
TAXO_LEVELS = ["domain", "phylum", "group", "order", "family", "specie"]
NB_LEVELS = 6
# nomenclature: {domain}_{phylum}_{group}_{order}_{family}_{specie}_{id}.fna  after merge_groupe_data.py
pattern_parts = [f"(?P<{level}>\\w+?)" for level in TAXO_LEVELS]
pattern_filename = "^" + "_".join(pattern_parts) + "(?:_(?P<id>\\d+))?\\.fna$"

def build_df_bacteria(datadir):
    data = []
    columns = TAXO_LEVELS + ["id"]

    for root, _, filenames in os.walk(datadir):
        for filename in filenames:
            path = os.path.join(root, filename)
            filename = os.path.basename(path).replace(".fna", "")
            parts = filename.split("_")
            if len(parts) == len(columns)==6:  # Vérification que le nombre de parties correspond
                data.append(parts)
            else:
                print(f"Format inattendu: {filename}")

    df_bacteria = pd.DataFrame(data, columns=columns)
    df_bacteria["id"] = df_bacteria["id"].astype(int)
    df_bacteria = df_bacteria.sort_values(by=TAXO_LEVELS[1:]+["id"], ascending=True).reset_index(drop=True)

    return df_bacteria


class BacteriaDataset:
    def __init__(self, datadir: str):

        if not os.path.isdir(datadir):
            raise ValueError(f"Le chemin {datadir} n'est pas un répertoire valide.")

        self.datadir = datadir

    def index_to_csv(self, csv_file):
        if os.path.exists(csv_file):
            print("RELOAD bacteria csv Index to ",csv_file)
            df_bacteria = pd.read_csv(csv_file, sep= ";")
            df_bacteria["id"] = df_bacteria["id"].astype(int)

        else:
            print("BUILD bacteriaReload csv Index to ", csv_file)
            df_bacteria = build_df_bacteria(self.datadir)
            df_bacteria.to_csv(csv_file, sep=";", index=False)

        return df_bacteria


if __name__ == '__main__':
    datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"
    dataset = BacteriaDataset(datadir=datadir)
    csv_file = "/home/hcourtei/Projects/MicroTaxo/codes/wisp/wisp_light/dataset/bacteria_index_csv"
    df_bacteria = dataset.index_to_csv(csv_file)
    # print(df_bacteria.to_markdown())
