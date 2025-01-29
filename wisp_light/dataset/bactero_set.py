import os

TAXO_LEVELS = ["domain", "phylum", "group", "order", "family", "specie"]

NB_LEVELS = 6

class BacteriaDataset:
    def __init__(self, datadir: str):

        if not os.path.isdir(datadir):
            raise ValueError(f"Le chemin {datadir} n'est pas un répertoire valide.")
        self.datadir = datadir

    def scan_all_files(self) :
        """
        Parcourt le répertoire de base et retourne la liste des fichiers.

        :return: Liste des chemins relatifs des fichiers.
        """
        files = []
        for root, _, filenames in os.walk(self.datadir):
            for filename in filenames:
                files.append(os.path.join(root, filename))

        return files

if __name__ == '__main__':
    datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_with_taxo"
    dataset = BacteriaDataset(datadir=datadir)
    files = dataset.scan_all_files()

    files_fna = [f  for f in files if f.endswith(".fna")]


    # tous les fichiers , sous forme de liste
    # filtrage fna + écriture autre dans non_fna.txt

    # merge des groupes ! atention group_0/nomA.fna  group_1/nomA.fna pour la sélection


    print(*files, sep='\n')