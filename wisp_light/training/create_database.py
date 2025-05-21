import json
import sys
from itertools import islice
from collections import Counter
from json import dumps
from functools import partial
from itertools import product
from dataclasses import dataclass
from Bio import SeqIO
from treelib import Tree
from concurrent.futures import ProcessPoolExecutor
from treelib.exceptions import DuplicatedNodeIdError
from tqdm import tqdm
from collections import defaultdict

from wisp_light.dataset.refSeqDataset import TAXO_LEVELS
from wisp_light.training.utils import log_resource_usage


def load_phylo_tree(database_json: str) :
    """
    Charge un arbre phylogénétique à partir d'une base de données JSON.

    Parameters
    ----------
    database_json : str
        Chemin vers le fichier JSON de la base de données.

    Returns
    -------
    Tree
        Arbre phylogénétique construit à partir de la base.
    int
        Nombre d'entrées de génomes indexées dans la base.
    """
    with open(database_json, 'r', encoding='utf-8') as f:
        db_data = json.load(f)  # Charge le fichier JSON

    # Initialisation de l'arbre
    phylo_tree = Tree()
    phylo_tree.create_node('Root', 'root_root', data=Taxonomy(0, 'Root', 'Root', None, None))

    json_datas = []  # Stocke les taxonomies pour recréer mapping_sp

    # Ajouter les nœuds en utilisant les taxonomies enregistrées
    for entry in db_data["datas"]:
        taxonomy_info = {level: entry.get(level, "Unknown") for level in TAXO_LEVELS}
        json_datas.append(taxonomy_info)

        # Construire la hiérarchie
        parent_id = "root_root"
        for level, name in taxonomy_info.items():
            node_id = f"{name.lower()}_{level}"
            if not phylo_tree.contains(node_id):  # Éviter les doublons
                phylo_tree.create_node(name, node_id, parent=parent_id, data=Taxonomy(None, level, name, None, None))
            parent_id = node_id  # Le parent du prochain niveau est le niveau actuel

    # Appliquer les codes aux nœuds selon les mappings
    taxa_codes = db_data["mappings"]
    for level in TAXO_LEVELS:
        if level in taxa_codes:
            for name, code in taxa_codes[level].items():
                node_id = f"{name.lower()}_{level}"
                if phylo_tree.contains(node_id):
                    phylo_tree.get_node(node_id).data.code = code

    # tree_output = phylo_tree.show(line_type="ascii", stdout=False)
    # logger.debug(tree_output)
    nb_genome_indexed = len(db_data['datas'])
    return phylo_tree, nb_genome_indexed

def build_database(train_dataset: list[str], params: dict, database_json: str , logger, max_workers) ->  Tree:
    """
    Construit une base de données phylogénétique et la sauvegarde au format JSON.

    Parameters
    ----------
    train_dataset : list of str
        Liste des tuples contenant le chemin vers le génome et sa taxonomie.
    params : dict
        Dictionnaire de paramètres (ksize, pattern, read_size, etc.).
    database_json : str
        Chemin du fichier JSON de sortie.
    logger : logging.Logger
        Objet de journalisation pour les messages d'information/débogage.
    max_workers : int
        Nombre de processus utilisés pour le calcul parallèle.

    Returns
    -------
    Tree
        Arbre phylogénétique résultant.
    """

    my_encoder: dict = encoder(ksize=params['ksize'])
    num_features = len(set(my_encoder.values()))

    # creating phylogenetic tree
    phylo_tree: Tree = Tree()
    phylo_tree.create_node('Root', 'root_root', data=Taxonomy(0, 'Root', 'Root', None, None))

    # Writing the database
    json_datas = list()
    count_kmers_partial = partial(counter_kmer, pattern=params['pattern'])   # Multiprocessing pour le comptage des kmers

    with open(database_json, 'w', encoding='utf-8') as jdb:
        jdb.write("{\n")
        jdb.write("\"datas\":[")
        # iterating over input genomes

        total_dna_length = 0
        total_nb_count_win = 0
        for id_genome, sample in tqdm(enumerate(train_dataset), disable=not sys.stdout.isatty()):
            genome_path, taxonomy = sample
            dna_sequence = read_genome(genome_path)

            # Splitting of reads
            total_dna_length += len(dna_sequence)
            # if len(dna_sequence) >= params['read_size']:
            all_reads = splitting(dna_sequence, params['read_size'], params['max_sampling'], shift_ratio=params['shift_ratio'])
            # l = list(all_reads)
            # Counting kmers inside each read
            with ProcessPoolExecutor(max_workers= max_workers) as executor:
                counters = list(executor.map(count_kmers_partial, all_reads))
                if id_genome== 10 :
                    log_resource_usage(logger, "build_database (Pool)")

            # counters: list[Counter] = [counter_kmer(read,params['ksize'],params['pattern']) for read in all_reads]
            total_nb_count_win += len(counters)
            del all_reads
            # Encoding reads for XGBoost
            encoded: list[dict] = [{my_encoder[k]:v for k, v in cts.items()} for cts in counters]
            del counters
            # Dumping in output file
            phylo_tree = taxonomy_information(taxonomy, phylo_tree)

            logger.debug("-"*60)
            logger.debug(f"\nTaxo {taxonomy}")
            tree_output = phylo_tree.show(line_type="ascii", stdout=False)
            logger.debug(tree_output)

            dict_write = {**taxonomy, 'datas': encoded}
            if id_genome == 0:
                jdb.write(dumps(dict_write))
            else :
                jdb.write(','+dumps(dict_write))
            del encoded

            json_datas.append({**taxonomy})  # if 'taxonomy' in locals(): #
            if 'taxonomy' not in locals():
                logger.error(f"ERROR NOT IN LOCALS {taxonomy}")
            del dna_sequence

        taxa_codes = mapping_sp(json_datas)
        jdb.write("],\"mappings\": ")         # Writing taxonomy to file
        jdb.write(dumps(taxa_codes))
        jdb.write("\n}")

        for i, level in enumerate(['root'] + TAXO_LEVELS):
            # list_node_depth_i = list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))
            node_iterator = (node for node in phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))
            logger.debug(f"Depth {i} level  {level}")
            for node in node_iterator:
                if node.data.code is None:
                    new_tag = taxa_codes[level][node.tag]
                    logger.debug(f"for NODE tag {node.tag} get node.data.code {new_tag}")
                    node.data.code = new_tag
        avg_dna_length = int(total_dna_length / len(train_dataset))
        avg_nb_count_win = int(total_nb_count_win / len(train_dataset))
        logger.info(f"avg_nb_count_win {avg_nb_count_win} avg_dna_length {avg_dna_length:,d}".replace(",", " "))
        log_resource_usage(logger, "build_database (end)")
    return phylo_tree


def mapping_sp(datas: list[dict]) -> dict:
    """
    Génère des codes de taxonomie pour chaque niveau à partir d’un ensemble de taxonomies.

    Parameters
    ----------
    datas : list of dict
        Liste de dictionnaires représentant les taxonomies pour chaque génome.

    Returns
    -------
    dict
        Dictionnaire des mappings : {niveau : {nom : code}}.
    """
    taxa_codes = defaultdict(lambda: {'number_taxa': 0})
    for sample in datas:
        for key, value in sample.items():
            if key in TAXO_LEVELS and value not in taxa_codes[key]:
                taxa_codes[key][value] = taxa_codes[key]['number_taxa']
                taxa_codes[key]['number_taxa'] += 1

    return taxa_codes


@dataclass
class Taxonomy:
    "Modélise un niveau taxonomique"
    code: int | None
    level: str
    name: str
    model_path: str | None
    config_path: str | None

def check_parameters(params: dict) :
    """
    Vérifie la validité des paramètres fournis.

    Parameters
    ----------
    params : dict
        Dictionnaire des paramètres à vérifier.

    Raises
    ------
    RuntimeError
        Si la somme des valeurs du motif n'est pas égale à `ksize`.
    KeyError
        Si la clé 'threshold' est absente des paramètres.
    """

    if not all([ sum(params['pattern']) == params['ksize'],] ):
        raise RuntimeError("Incorrect parameter file")
    if "threshold"  not in params:
        raise KeyError("Invalid parameter file, must contain a read acceptance threshold value "
            " between 0.01 (1% identity) and 1.0 (100% identity).")
    if params["threshold"] > 1.0:
        params["threshold"] = 1.0
    elif  params["threshold"] < 0.01:
        params["threshold"] = 0.01

    # return  # verifies that the pattern length respects ksize



def encoder(ksize: int) -> dict:
    """
    Génère un encodage k-mer mappant les chaînes d'ADN à des codes entiers.

    Parameters
    ----------
    ksize : int
       Longueur du k-mer.

    Returns
    -------
    dict
       Dictionnaire mappant les chaînes k-mer aux entiers.
    """

    res = {code: encode_kmer(code) for code in map(''.join, product('ATCG', repeat=ksize))}
    return res


def encode_kmer(kmer: str) -> int:
    """
    Encode une chaîne k-mer en un entier en base 4.

    Parameters
    ----------
    kmer : str
        K-mer ADN composé de A, T, C, G.

    Returns
    -------
    int
        Encodage entier en base 4 du k-mer.
    """
    mapper: dict = {'A': "0",'C': "1", 'G': "2", 'T': "3",}
    return int(''.join([mapper[k] for k in kmer]))


def taxonomy_information(taxo_dict:dict, tree_struct: Tree) ->  Tree:
    """
    Insère les nœuds de taxonomie dans une structure d'arbre.

    Parameters
    ----------
    taxo_dict : dict
        Informations taxonomiques {niveau : nom}.
    tree_struct : Tree
        Arbre existant à mettre à jour.

    Returns
    -------
    Tree
        Arbre mis à jour avec les nœuds de taxonomie.
    """
    # genome_path = "Bacteria_Pseudomonadati_Bacteroidota_Flavobacteriales_Elizabethkingia_meningoseptica.fna"
    # taxa_family_old = ['root'] + TAXO_LEVELS[1:]
    assert list(taxo_dict.keys()) == TAXO_LEVELS
    taxa_family = ['root'] + list(taxo_dict.keys())
    # only_name_file = os.path.splitext(os.path.basename(genome_path))[0]
    # example 'Bacteria_Campylobacterota_Epsilonproteobacteria_Campylobacterales_Campylobacter_hyointestinalis'
    # taxo_info: list = only_name_file.split('_')[:5] # only 5 first
    # ['Bacteria', 'Campylobacterota', 'Epsilonproteobacteria', 'Campylobacterales', 'Campylobacter']
    # parents_old = ['Root'] + taxo_info
    parents = ['Root'] + list(taxo_dict.values())

    for i, x in enumerate(parents):
        # y = x
        if x != 'Root':
            idx = f"{x.lower()}_{taxa_family[i]}"
            parent = f"{parents[i - 1].lower()}_{taxa_family[i - 1]}"
            level = taxa_family[1:][i - 1]
            name = x
            data = Taxonomy(None, level, name, None, None)
            try:
                tree_struct.create_node(x, idx, parent=parent, data=data)
            except DuplicatedNodeIdError as e :
                # logger.error(f"Error: Duplicate node with  {os.path.basename(genome_path)}.  {e}")
                pass

    return tree_struct # dict(zip(TAXO_LEVELS, taxo_info)),

ALL_COMPLEMENTS =  { 'A': 'T',
                     'T': 'A',
                     'C': 'G',
                     'G': 'C',
                     'U': 'A',
                     'R': 'Y',
                     'Y': 'R',
                     'K': 'M',
                     'M': 'K',
                     'S': 'W',
                     'W': 'S',
                     'B': 'V',
                     'V': 'B',
                     'D': 'H',
                     'H': 'D',
                     'N': 'N'}
NUCLEOTIDE_MAP = {
    'A': ['A'], 'T': ['T'], 'C': ['C'], 'G': ['G'],
    'U': ['T'],  # Uracile → Thymine
    'R': ['G', 'A'], 'Y': ['C', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'S': ['G', 'C'], 'W': ['A', 'T'], 'B': ['G', 'T', 'C'],
    'D': ['G', 'T', 'A'], 'H': ['A', 'T', 'C'], 'V': ['G', 'A', 'C']
}

# Par défaut, un caractère inconnu peut être A, T, C ou G
DEFAULT_NUC = ['A', 'T', 'C', 'G']



def pattern_filter(substring: str, pattern: list[int]) -> str:
    """
    Applique un filtrage basé sur un motif à une sous-chaîne.

    Parameters
    ----------
    substring : str
        Sous-chaîne à filtrer.
    pattern : list of int
        Liste indiquant la répétition des caractères.

    Returns
    -------
    str
        Chaîne filtrée selon le motif.
    """
    pattern_filtered = ''.join([char * pattern[i] for i, char in enumerate(substring)])
    return pattern_filtered


def revcomp(seq: str) -> str:
    """
    Calcule le complément inverse d'une séquence ADN.

    Parameters
    ----------
    seq : str
        Séquence ADN.

    Returns
    -------
    str
        Complément inverse de la séquence d'entrée.
    """
    # reverse_complement = ''.join([compl[s] for s in string][::-1]) original
    # reverse_complement = ''.join(map(ALL_COMPLEMENTS.get, reversed(string)))
    reverse_complement =''.join(ALL_COMPLEMENTS[c] for c in reversed(seq))
    return reverse_complement #seq.translate(COMP_TABLE)[::-1]

#
def process_ambiguous(counts):
    """
    Développe les k-mers contenant des bases ambiguës en toutes leurs combinaisons possibles.

    Parameters
    ----------
    counts : Counter or dict
        Dictionnaire de comptage de k-mers, pouvant contenir des bases ambiguës.

    Returns
    -------
    Counter
        Comptage mis à jour avec les ambiguïtés résolues.
    """
    counts = defaultdict(int, counts)  # Convertir counts en defaultdict(int)
    modifications = []  # Stocker les nouvelles valeurs
    to_remove = []  # Stocker les clés à supprimer

    for key, count in counts.items():
        list_of_keys = ['']
        for x in key:
            nuct = NUCLEOTIDE_MAP.get(x, DEFAULT_NUC)  # Liste des nucléotides possibles
            list_of_keys = [new_key + n for new_key in list_of_keys for n in nuct]

            # Normalisation : Répartir le compte entre toutes les combinaisons
        kmer_number = count // len(list_of_keys)
        modifications.extend((prob_key, kmer_number) for prob_key in list_of_keys)
        to_remove.append(key)  # On supprime l'ancienne clé après l’itération

    # Supprimer les anciennes clés en toute sécurité
    for key in to_remove:
        counts.pop(key, None)  # pop() évite les erreurs si la clé est absente

    # Appliquer les nouvelles valeurs
    for key, value in modifications:
        counts[key] += value  # defaultdict évite .get(), code plus propre et rapide

    return counts

#
# Fonction de lecture
def read_genome(genome_path):
    """
    Lit une séquence génomique depuis un fichier FASTA.

    Parameters
    ----------
    genome_path : str
        Chemin vers le fichier génomique.

    Returns
    -------
    str
        Séquence ADN complète en une seule chaîne.
    """

    with open(genome_path, 'r', encoding='utf-8') as freader:
        # genome_data = (str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta'))
        return ''.join(str(fasta.seq).upper() for fasta in SeqIO.parse(freader, 'fasta'))


def splitting(seq: str, read_size: int, max_sampling = None, shift_ratio = None) :
    """
    Découpe une séquence génomique en sous-fenêtres.

    Parameters
    ----------
    seq : str
        Séquence ADN à découper.
    read_size : int
        Longueur de chaque sous-fenêtres.
    max_sampling : int, optional
        Nombre maximum de sous-fenêtres.
    shift_ratio : float, optional
        Ratio de recouvrement entre les sous-fenêtres.

    Returns
    -------
    iterator
        Itérateur sur les sous-fenêtres.

    Raises
    ------
    ValueError
        Si la sous-fenêtre est trop courte ou si les paramètres sont incohérents.
    """
    if len(seq) < read_size:
        raise ValueError("Read is too short.")
    if shift_ratio is not None and max_sampling is not None:
        raise ValueError("Provide either shift_ratio or max_sampling, not both.")

    shift = int(shift_ratio * read_size) if shift_ratio else int((len(seq)-read_size)/max_sampling)
    if shift <= 0 :
        raise ValueError(f"Calculated shift {shift} must be positive.")
    # for i in range(nb_win):
    #     yield seq[shift*i:shift*i+read_size]
    return islice((seq[i:i + read_size] for i in range(0, len(seq) - read_size + 1, shift)), max_sampling)

def counter_kmer(read: str, pattern: list[int]) -> Counter:
    """
    Compte les k-mers dans une lecture, y compris les compléments inverses et les ambiguïtés.

    Parameters
    ----------
    read : str
        Lecture ADN.
    pattern : list of int
        Motif utilisé pour filtrer ou transformer les k-mers.

    Returns
    -------
    Counter
        Dictionnaire des k-mers comptés après traitement.
    """
    kmer_size = len(pattern)

    all_kmers = (read[i:i+kmer_size] for i in range(len(read)-kmer_size+1))
    # counts = Counter(map(lambda i: read[i:i + kmer_size], range(len(read) - kmer_size + 1))) ## --
    counts = Counter(all_kmers)
    counts +=  Counter({revcomp(k): v for k, v in counts.items()}) #comptage avec reverse-complements.


    if not all(pattern):
        # All positions in pattern should not be kept, we apply filter
        counts = Counter({pattern_filter(k, pattern): v for k, v in counts.items()})

    # suppression des homopolymères (AAAA, TTTT, etc.)
    for base in "ATCG":
        counts.pop(base * kmer_size, None)

        # We treat cases where sequence alphabet is not ATCG, Gérer les nucléotides ambigus (RYKM...U).
    # counts = dict(counts)  # Correction ici
    counts_purged = process_ambiguous(counts)

    return Counter(counts_purged)