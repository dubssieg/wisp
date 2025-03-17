import logging
import os
import sys
import pandas as pd
from sklearn.model_selection import train_test_split
from tqdm import tqdm
sys.path.append('../../..')

TAXO_LEVELS =  ['phylum', 'class', 'order', 'family']

class RefSeqDataset:
    def __init__(self, index_csv, datadir, logger=None):
        self.index_with_label = pd.read_csv(index_csv, sep='\t')
        self.datadir = datadir
        self.logger = logger
        self.pairing_label_to_file()

    def pairing_label_to_file(self):
        all_files = os.listdir(self.datadir)

        gcf_ids_from_labels = self.index_with_label['Assembly Accession']
        file_map = {gcf_id: next((file for file in all_files if gcf_id in file), None) for gcf_id in gcf_ids_from_labels}

        self.index_with_label['file'] = self.index_with_label['Assembly Accession'].map(file_map)
        self.index_with_label = self.index_with_label.dropna(subset=['file']).reset_index(drop=True)
        self.logger.info(f"nb files in datadir: {len(all_files)} restrict to {len(self.index_with_label)} with labels in index ")

    def __len__(self):
        # print(f"Current dataset length: {len(self.index_with_label)}")
        return len(self.index_with_label)

    def __getitem__(self, idx):
        row = self.index_with_label.loc[idx]
        file_fna = os.path.join(self.datadir, row['file'])
        taxonomy_info = self.index_with_label.loc[idx,TAXO_LEVELS].to_dict()
        return file_fna, taxonomy_info

    def get_indices(self):
        """Retourne les indices du DataFrame index_with_label."""
        return self.index_with_label.index

    def split(self, test_size=0.2, random_state=None, family_strat=False):
        """Effectue un split aléatoire et retourne deux instances de RefSeqDataset."""
        strat = None if not family_strat  else self.index_with_label['family']
        train_df, test_df = train_test_split(self.index_with_label, test_size=test_size,
                                             random_state=random_state, stratify=strat)
        train_df = train_df.reset_index(drop=True)
        test_df = test_df.reset_index(drop=True)
        # Créer deux nouvelles instances de RefSeqDataset pour les ensembles d'entraînement et de test
        train_dataset = RefSeqDataset.from_dataframe(train_df, self.datadir)
        test_dataset = RefSeqDataset.from_dataframe(test_df, self.datadir)
        self.logger.info(f"Splitted dataset nb {len(self.index_with_label)} into train :{len(train_dataset)} val: {len(test_dataset)}")

        return train_dataset, test_dataset

    @classmethod
    def from_dataframe(cls, dataframe, datadir):
        """Créer une instance de RefSeqDataset à partir d'un DataFrame existant."""
        instance = cls.__new__(cls)
        instance.index_with_label = dataframe
        instance.datadir = datadir
        return instance

    def __iter__(self):
        """Permet l'itération sur les éléments de l'index_with_label."""
        for idx in range(len(self.index_with_label)):
            yield self[idx]





if __name__ == '__main__':
    from wisp.wisp_light.training.utils import setup_logger

    logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

    index_csv = 'complete_refseq_referent_genome_with_taxo.tsv'
    datadir = '/projects/microtaxo/data/refseq3'
    # datadir = '/home/hcourtei/Projects/MicroTaxo/codes/data/refseq/group_1'
    ds = RefSeqDataset(index_csv,datadir, logger)

    # print(ds.index_with_label['file'])
    x, y = ds[0]

    train_dataset, test_dataset = ds.split(test_size=0.1, random_state=42, family_strat=False)
    # print("lenght train_dataset", len(train_dataset))
    # print("index", train_dataset.index_with_label.index)


    for id_genome, sample in enumerate(train_dataset):
        if id_genome >4:
            continue
        genome, taxo_dict = sample
        print(f" {id_genome}, genome {genome}")
        print("taxo_dict", taxo_dict)

    # split sur les indices
    # indices = ds.get_indices()
    # train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
    # # Accéder aux ensembles d'entraînement et de test via les indices
    # train_set = ds.index_with_label.loc[train_indices]
    # test_set = ds.index_with_label.loc[test_indices]
    #
    # # Afficher les tailles des ensembles
    # print(f"Training set size: {len(train_set)}")
    # print(f"Test set size: {len(test_set)}")


    # Créer un KFold avec 5 splits
    # from sklearn.model_selection import KFold
    # kfold = KFold(n_splits=5, shuffle=True, random_state=42)
    #
    # # Itérer sur chaque fold
    # for fold, (train_indices, test_indices) in enumerate(kfold.split(indices)):
    #     print(f"Fold {fold + 1}")
    #
    #     # Sélectionner les ensembles d'entraînement et de test à partir des indices
    #     train_set = ds.index_with_label.loc[train_indices]
    #     test_set = ds.index_with_label.loc[test_indices]
    #
    #     # Afficher les tailles des ensembles pour ce fold
    #     print(f"Training set size: {len(train_set)}")
    #     print(f"Test set size: {len(test_set)}")