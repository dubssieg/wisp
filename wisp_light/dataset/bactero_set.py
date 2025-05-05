import logging
import os
import re
import pandas as pd
from sklearn.model_selection import train_test_split

from wisp_light.training.utils import setup_logger

# Constantes de niveaux taxonomiques
TAXO_LEVELS = ["domain", "phylum", "group", "order", "family", "specie"]
NB_LEVELS = 6

# Expression régulière pour parser les fichiers de séquences .fna
# Format attendu : domain_phylum_group_order_family_specie_id.fna
pattern_parts = [f"(?P<{level}>\\w+?)" for level in TAXO_LEVELS]
pattern_filename = "^" + "_".join(pattern_parts) + "_(?P<id>\\d+)\\.fna$"


class BacteriaDataset:
    """
    Classe représentant un jeu de données de bactéries structuré par niveaux taxonomiques.

    Parameters
    ----------
    datadir : str
        Chemin vers le répertoire contenant les fichiers `.fna`.
    logger : logging.Logger, optional
        Logger pour la journalisation des événements. Si None, un logger par défaut est créé.

    Attributes
    ----------
    datadir : str
        Répertoire contenant les fichiers d'entrée.
    logger : logging.Logger
        Logger configuré.
    csv_file : str
        Chemin vers le fichier CSV d’index des bactéries.
    df_bacteria : pd.DataFrame
        DataFrame contenant les métadonnées des fichiers bactériens.
    df_selected : pd.DataFrame or None
        DataFrame contenant les données filtrées selon les critères choisis.
    """

    def __init__(self, datadir: str, logger=None):
        if logger is None:
            logger = setup_logger(os.path.basename(__file__), level=logging.INFO, log_file=None)

        if not os.path.isdir(datadir):
            raise ValueError(f"Le chemin {datadir} n'est pas un répertoire valide.")

        self.datadir = datadir
        self.logger = logger
        self.csv_file = os.path.join(os.path.dirname(datadir), 'bacteria_index.csv')
        self.df_bacteria = self.build()
        self.df_selected = None

    def build(self):
        """
        Construit ou recharge l'index des bactéries à partir des fichiers `.fna`.

        Returns
        -------
        pd.DataFrame
            DataFrame contenant les informations taxonomiques extraites des noms de fichiers.
        """
        if os.path.exists(self.csv_file):
            self.logger.info(f"RELOAD bacteria csv Index from {self.csv_file}")
            df_bacteria = pd.read_csv(self.csv_file, sep=";")
            df_bacteria["id"] = df_bacteria["id"].astype(int)
        else:
            self.logger.info(f"BUILD bacteria csv Index to {self.csv_file}")
            df_bacteria = build_df_bacteria(self.datadir, self.logger)
            df_bacteria.to_csv(self.csv_file, sep=";", index=False)

        return df_bacteria

    def filter_family_by_min_species(self, min_family_threshold=5, max_family_repr: int | str = "inf"):
        """
        Filtre les familles ayant au moins un certain nombre d'espèces représentées.

        Parameters
        ----------
        min_family_threshold : int
            Nombre minimum d'espèces différentes par famille à retenir.
        max_family_repr : int or str
            Nombre maximum de représentants à conserver par famille. Peut être "inf" pour ignorer.
        """
        df_representants = self.df_bacteria[self.df_bacteria["id"] == 0].reset_index(drop=True)

        result = df_representants.groupby(['domain', 'phylum', 'group', 'order', 'family']) \
            ['specie'].nunique().reset_index(name='specie_count')

        result_sorted = result.sort_values(by='specie_count', ascending=False).reset_index(drop=True)

        filtered_result = result[result['specie_count'] >= min_family_threshold]

        self.df_selected = df_representants.merge(
            filtered_result[['domain', 'phylum', 'order', 'family']],
            on=['domain', 'phylum', 'order', 'family'],
            how='inner'
        )

        self.logger.info(f"Selected .fna for families with >= {min_family_threshold} species.")
        self.logger.info(f"initial genomes: {len(self.df_bacteria)}, "
                         f"representants: {len(df_representants)}, "
                         f"final filter: {len(self.df_selected)}")

        if float(max_family_repr) < float("inf"):
            assert int(max_family_repr) > min_family_threshold
            self.df_selected = self.df_selected.groupby("family").head(int(max_family_repr))
            self.logger.info(f"Selected .fna with max {max_family_repr} representatives per family")

    def train_test_split(self, test_size=0.2, random_state=42):
        """
        Effectue un découpage train/test tout en conservant la distribution par famille.

        Parameters
        ----------
        test_size : float
            Proportion du jeu de test (ex: 0.2 pour 20%).
        random_state : int
            Graine pour la reproductibilité.

        Returns
        -------
        train_files_list : list of str
            Chemins complets vers les fichiers d'entraînement.
        val_files_list : list of str
            Chemins complets vers les fichiers de validation.
        """
        if hasattr(self, 'df_selected'):
            self.df_selected['filename'] = self.df_selected.apply(
                lambda row: f"{row['domain']}_{row['phylum']}_{row['group']}_{row['order']}_"
                            f"{row['family']}_{row['specie']}_{row['id']}.fna", axis=1
            )

            df = self.df_selected[['filename', 'family']]

            # Split par famille
            train_files, val_files = zip(*[
                train_test_split(group['filename'], test_size=test_size, random_state=random_state)
                for _, group in df.groupby('family')
            ])

            # Reconstituer les DataFrames
            df_train = df[df['filename'].isin([f for files in train_files for f in files])]
            df_val = df[df['filename'].isin([f for files in val_files for f in files])]

            train_files_list = [os.path.join(self.datadir, fname) for fname in df_train['filename'].tolist()]
            val_files_list = [os.path.join(self.datadir, fname) for fname in df_val['filename'].tolist()]

            self.logger.info(f"Splitted dataset into train: {len(train_files_list)} and val: {len(val_files_list)}")

        return train_files_list, val_files_list


def build_df_bacteria(datadir, logger):
    """
    Construit un DataFrame en extrayant les informations taxonomiques depuis les noms de fichiers `.fna`.

    Parameters
    ----------
    datadir : str
        Répertoire contenant les fichiers `.fna`.
    logger : logging.Logger
        Logger pour journaliser les erreurs ou avertissements.

    Returns
    -------
    pd.DataFrame
        DataFrame avec colonnes taxonomiques et identifiants.
    """
    data = []
    columns = TAXO_LEVELS + ["id"]
    df_bacteria = pd.DataFrame(data, columns=columns)

    for root, _, filenames in os.walk(datadir):
        for filename in filenames:
            match = re.match(pattern_filename, filename)
            if match:
                result = match.groupdict()
                df_bacteria.loc[len(df_bacteria)] = result
            else:
                logger.info(f"WARNING: filename does not match pattern: {filename}")

    df_bacteria["id"] = df_bacteria["id"].astype(int)
    df_bacteria = df_bacteria.sort_values(by=TAXO_LEVELS[1:] + ["id"], ascending=True).reset_index(drop=True)

    return df_bacteria


if __name__ == '__main__':
    datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo_merged"

    dataset = BacteriaDataset(datadir)
    dataset.filter_family_by_min_species(min_family_threshold=5, max_family_repr='inf')
    print(dataset.df_selected.to_markdown())
