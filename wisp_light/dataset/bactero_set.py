import logging
import os
import re
import pandas as pd
from sklearn.model_selection import train_test_split
import sys

sys.path.append('..')
from training.utils import setup_logger

TAXO_LEVELS = ["domain", "phylum", "group", "order", "family", "specie"]
NB_LEVELS = 6
# nomenclature: {domain}_{phylum}_{group}_{order}_{family}_{specie}_{id}.fna  after merge_groupe_data.py
pattern_parts = [f"(?P<{level}>\\w+?)" for level in TAXO_LEVELS]
pattern_filename = "^" + "_".join(pattern_parts) + "_(?P<id>\\d+)\\.fna$" # A_B_C_D_E_F_456.fna  un nombre obligatoire
# regex = re.compile(pattern_filename)

class BacteriaDataset:
    def __init__(self, datadir: str, logger=None):

        if logger is  None:
            logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

        if not os.path.isdir(datadir):
            raise ValueError(f"Le chemin {datadir} n'est pas un répertoire valide.")

        self.datadir = datadir
        self.logger = logger
        self.csv_file = os.path.join(os.path.dirname(datadir), 'bacteria_index.csv') # "/home/hcourtei/Projects/MicroTaxo/codes/wisp/wisp_light/dataset/"
        self.df_bacteria = self.build()
        self.df_selected = None

    def build(self):
        if os.path.exists(self.csv_file):
            self.logger.info(f"RELOAD bacteria csv Index to {self.csv_file}", )
            df_bacteria = pd.read_csv(self.csv_file, sep= ";")
            df_bacteria["id"] = df_bacteria["id"].astype(int)

        else:
            self.logger.info("BUILD bacteriaReload csv Index to ", self.csv_file)
            df_bacteria = build_df_bacteria(self.datadir, self.logger)
            df_bacteria.to_csv(self.csv_file, sep=";", index=False)

        return df_bacteria

    def filter_family_by_min_species(self, min_family_threshold = 5, max_family_repr: int |str  = "inf"):

        df_representants = self.df_bacteria[self.df_bacteria["id"] == 0].reset_index(drop=True)
        result = df_representants.groupby(['domain', 'phylum', 'group', 'order', 'family']) \
            ['specie'].nunique().reset_index(name='specie_count')

        result_sorted = result.sort_values(by='specie_count', ascending=False).reset_index(drop=True)
        filtered_result = result[(result['specie_count'] >= min_family_threshold)]# &
        #                          (result['specie_count'] < max_family_threshold)]
        self.df_selected = df_representants.merge(filtered_result[['domain', 'phylum', 'order', 'family']],
                                             on=['domain', 'phylum', 'order', 'family'],
                                             how='inner')
        self.logger.info(f"Selected .fna for family with nb differents species >= {min_family_threshold}")
        self.logger.info(f"initial genome: {len(self.df_bacteria)} representants {len(df_representants)} to final filter {len(self.df_selected)}")
        if float(max_family_repr) < float("inf"):
            assert max_family_repr > min_family_threshold
            self.df_selected = self.df_selected.groupby("family").head(max_family_repr)
            self.logger.info(f"Selected .fna for family with max representant {max_family_repr}")


    def train_test_split(self, test_size=0.2, random_state=42):
        if hasattr(self, 'df_selected'):

            self.df_selected['filename'] = self.df_selected.apply(
                lambda
                    row: f"{row['domain']}_{row['phylum']}_{row['group']}_{row['order']}_{row['family']}_{row['specie']}_{row['id']}.fna",
                axis=1
            )
            df = self.df_selected[['filename', 'family']]

            # Split 80%/20% pour chaque famille
            train_files, val_files = zip(*[
                train_test_split(group['filename'], test_size=test_size, random_state=random_state)
                for _, group in df.groupby('family')
            ])

            # Rassembler les résultats dans les DataFrames
            df_train = df[df['filename'].isin([file for files in train_files for file in files])]
            df_val = df[df['filename'].isin([file for files in val_files for file in files])]

            train_files_list = df_train['filename'].tolist()
            val_files_list = df_val['filename'].tolist()
            train_files_list = [os.path.join(self.datadir, filename) for filename in train_files_list]
            val_files_list = [os.path.join(self.datadir, filename) for filename in val_files_list]
            self.logger.info(f"Splitted dataset into train :{len(train_files_list)} val: {len(val_files_list)}")

        return train_files_list, val_files_list




def build_df_bacteria(datadir, logger):
    data = []
    columns = TAXO_LEVELS + ["id"]
    df_bacteria = pd.DataFrame(data, columns=columns)

    for root, _ , filenames in os.walk(datadir):
        for filename in filenames:
            match = re.match(pattern_filename, filename)
            if match:
                result = match.groupdict()
                df_bacteria.loc[len(df_bacteria)] = result
            else:  # if no match pattern
                logger.info(f"WARNING no pattern match for {filename}")

    df_bacteria["id"] = df_bacteria["id"].astype(int)
    df_bacteria = df_bacteria.sort_values(by=TAXO_LEVELS[1:]+["id"], ascending=True).reset_index(drop=True)

    return df_bacteria




if __name__ == '__main__':

    datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"

    dataset = BacteriaDataset(datadir)
    # print(dataset.df_bacteria.to_markdown())
    dataset.filter_family_by_min_species(min_family_threshold = 5, max_family_repr='inf')
    print(dataset.df_selected.to_markdown())

    # print(*file_names, sep='\n')
    # train_files_list, val_files_list =  dataset.train_test_split(test_size=0.2, random_state=42)

    # # Afficher les résultats
    # print("Ensemble d'entraînement:")
    # print(*train_files_list, sep='\n')
    # print("\nEnsemble de validation:")
    # print(*val_files_list, sep='\n')
