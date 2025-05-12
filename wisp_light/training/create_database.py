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
    Loads a phylogenetic tree from a JSON database.

    Parameters
    ----------
    database_json : str
        Path to the JSON database file.

    Returns
    -------
    tuple
        Tree : `Tree`
            A phylogenetic tree constructed from the database.
        int
            Number of genome entries indexed in the database.
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
    Builds a phylogenetic database and saves it as JSON.

    Parameters
    ----------
    train_dataset : list of str
        List of tuples containing genome path and taxonomy.
    params : dict
        Dictionary of parameters including ksize, pattern, read_size, etc.
    database_json : str
        Path to the output JSON file.
    logger : logging.Logger
        Logger for debug and info output.
    max_workers : int
        Number of workers for multiprocessing.

    Returns
    -------
    Tree
        The constructed phylogenetic tree.
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
    Generates taxonomy codes per level from a dataset.

    Parameters
    ----------
    datas : list of dict
        List of taxonomy dictionaries per genome.

    Returns
    -------
    dict
        Dictionary of mappings {level: {name: code}}.
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
    "Modelizes a taxa level"
    code: int | None
    level: str
    name: str
    model_path: str | None
    config_path: str | None

def check_parameters(params: dict) :
    """
    Checks the validity of input parameters.

    Parameters
    ----------
    params : dict
        Dictionary of parameters to validate.

    Raises
    ------
    RuntimeError
        If the sum of pattern values is not equal to ksize.
    KeyError
        If 'threshold' is not present in params.
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
    Generates a k-mer encoder mapping DNA strings to integer codes.

    Parameters
    ----------
    ksize : int
        Length of the k-mer.

    Returns
    -------
    dict
        Mapping of k-mer strings to integer codes.
    """

    res = {code: encode_kmer(code) for code in map(''.join, product('ATCG', repeat=ksize))}
    return res


def encode_kmer(kmer: str) -> int:
    """
    Encodes a k-mer string into base-4 integer representation.

    Parameters
    ----------
    kmer : str
        DNA k-mer composed of A, T, C, G.

    Returns
    -------
    int
        Base-4 encoded integer of the k-mer.
    """
    mapper: dict = {'A': "0",'C': "1", 'G': "2", 'T': "3",}
    return int(''.join([mapper[k] for k in kmer]))


def taxonomy_information(taxo_dict:dict, tree_struct: Tree) ->  Tree:
    """
    Inserts taxonomy nodes into a tree structure.

    Parameters
    ----------
    taxo_dict : dict
        Taxonomic information {level: name}.
    tree_struct : Tree
        Existing tree structure to update.

    Returns
    -------
    Tree
        Updated tree structure.
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
    Applies a pattern-based character filtering to a substring.

    Parameters
    ----------
    substring : str
        Substring to filter.
    pattern : list of int
        List indicating character repetition.

    Returns
    -------
    str
        Filtered string according to the pattern.
    """
    pattern_filtered = ''.join([char * pattern[i] for i, char in enumerate(substring)])
    return pattern_filtered


def revcomp(seq: str) -> str:
    """
    Computes the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Reverse complement of the input sequence.
    """
    # reverse_complement = ''.join([compl[s] for s in string][::-1]) original
    # reverse_complement = ''.join(map(ALL_COMPLEMENTS.get, reversed(string)))
    reverse_complement =''.join(ALL_COMPLEMENTS[c] for c in reversed(seq))
    return reverse_complement #seq.translate(COMP_TABLE)[::-1]

#
def process_ambiguous(counts):
    """
    Expands ambiguous nucleotide kmers into their possible combinations.

    Parameters
    ----------
    counts : Counter or dict
        K-mer count dictionary, possibly containing ambiguous bases.

    Returns
    -------
    Counter
        Updated k-mer counts with ambiguity resolved.
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
   Reads a genome sequence from a FASTA file.

   Parameters
   ----------
   genome_path : str
       Path to the genome file.

   Returns
   -------
   str
       Full DNA sequence as a single string.
   """

    with open(genome_path, 'r', encoding='utf-8') as freader:
        # genome_data = (str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta'))
        return ''.join(str(fasta.seq).upper() for fasta in SeqIO.parse(freader, 'fasta'))


def splitting(seq: str, read_size: int, max_sampling = None, shift_ratio = None) :
    """
    Splits a genome sequence into subreads.

    Parameters
    ----------
    seq : str
       DNA sequence to split.
    read_size : int
       Length of each subread.
    max_sampling : int, optional
       Maximum number of subreads.
    shift_ratio : float, optional
       Ratio of overlap between reads.

    Returns
    -------
    iterator
       Iterator over subreads.

    Raises
    ------
    ValueError
       If the read is too short or parameters are inconsistent.
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
    Counts kmers in a read, including reverse complements and ambiguous resolutions.

    Parameters
    ----------
    read : str
       DNA read.
    pattern : list of int
       Pattern used to filter or transform k-mers.

    Returns
    -------
    Counter
       Count of kmers after filtering and resolution.
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