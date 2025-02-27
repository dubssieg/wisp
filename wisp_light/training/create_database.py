"""Creates a json database"""
import json
import logging
import os
from collections import Counter
from json import dumps
from typing import Generator
from itertools import product
from dataclasses import dataclass
from Bio import SeqIO
from treelib import Tree
from treelib.exceptions import DuplicatedNodeIdError
from tqdm import tqdm
from utils import setup_logger
from collections import defaultdict

TAXO_LEVELS = ["domain", "phylum", "group", "order", "family"]




def load_phylo_tree(databse_json: str) :
    with open(databse_json, 'r', encoding='utf-8') as f:
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

def build_database(input_file_list: list[str], params: dict, database_json: str , logger) ->  Tree:
    """Builds a json file with taxa levels as dict information"""
    # creating encoder
    my_encoder: dict = encoder(ksize=params['ksize'])
    # creating phylogenetic tree
    phylo_tree: Tree = Tree()
    phylo_tree.create_node('Root', 'root_root', data=Taxonomy(0, 'Root', 'Root', None, None))

    # Writing the database
    json_datas = list()

    with open(database_json, 'w', encoding='utf-8') as jdb:
        jdb.write("{\n")
        jdb.write("\"datas\":[")
        # iterating over input genomes

        total_dna_length = 0
        total_nb_count_win = 0
        for id_genome, genome in (pbar:= tqdm(enumerate(input_file_list))):
            # pbar.set_description(f"Genome {path.basename(genome)}")
            with open(genome, 'r', encoding='utf-8') as freader:
                # genome_data: list = [str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta')]
                genome_data = (str(fasta.seq) for fasta in SeqIO.parse(freader, 'fasta'))  # Générateur
                dna_sequence = (''.join([seq for seq in genome_data])).upper() # Merging all seqs together
            # Splitting of reads
            total_dna_length += len(dna_sequence)
            if len(dna_sequence) >= params['read_size']:
                all_reads = splitting(dna_sequence, params['read_size'], params['max_sampling'], shift_ratio=params['shift_ratio'])
                # l = list(all_reads)
                # Counting kmers inside each read
                counters: list[Counter] = [counter_kmer(read,params['ksize'],params['pattern']) for read in all_reads]
                total_nb_count_win += len(counters)
                del all_reads
                # Encoding reads for XGBoost
                encoded: list[dict] = [{my_encoder[k]:v for k, v in cts.items()} for cts in counters]
                del counters
                # Dumping in output file
                taxonomy, phylo_tree = taxonomy_information(genome, phylo_tree)

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
            del genome_data

        taxa_codes = mapping_sp(json_datas)
        jdb.write("],\"mappings\": ")         # Writing taxonomy to file
        jdb.write(dumps(taxa_codes))
        jdb.write("\n}")
        # numérotation des noeuds de l'arbre
        # Depth 2 level  phylum
        # from NODE tag Campylobacterota   to  0
        # from NODE tag Spirochaetota   to  1
        # from NODE tag Pseudomonadati   to  2
        for i, level in enumerate(['root'] + TAXO_LEVELS):
            # list_node_depth_i = list(phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))
            node_iterator = (node for node in phylo_tree.filter_nodes(lambda x: phylo_tree.depth(x) == i))
            logger.debug(f"Depth {i} level  {level}")
            for node in node_iterator:
                if node.data.code is None:
                    new_tag = taxa_codes[level][node.tag]
                    logger.debug(f"for NODE tag {node.tag} get node.data.code {new_tag}")
                    node.data.code = new_tag
        avg_dna_length = int(total_dna_length / len(input_file_list))
        avg_nb_count_win = int(total_nb_count_win / len(input_file_list))
        logger.info(f"avg_nb_count_win {avg_nb_count_win} avg_dna_length {avg_dna_length:,d}".replace(",", " "))
    return phylo_tree


def mapping_sp(datas: list[dict]) -> dict:
    """Given a dataset, creates numbers for each level of classification
    Args:
        datas (list[dict]): a set of assigned reads
    Returns:
        dict: a grouped-by-level list of codes
    """
    taxa_codes = defaultdict(lambda: {'number_taxa': 0})
    for sample in datas:
        for key, value in sample.items():
            if key in TAXO_LEVELS and value not in taxa_codes[key]:
                taxa_codes[key][value] = taxa_codes[key]['number_taxa']
                taxa_codes[key]['number_taxa'] += 1


    # taxa_codes: dict = {taxon: {'number_taxa': 0} for taxon in TAXO_LEVELS}
    # for sample in datas:
    #     for key, value in sample.items():
    #         if key in TAXO_LEVELS:
    #             taxa_level = taxa_codes[key]
    #             if not value in taxa_level:
    #                 taxa_level[value] = taxa_level['number_taxa']
    #                 taxa_level['number_taxa'] += 1
    return taxa_codes


def splitting(seq: str, read_size: int, max_sampling = None, shift_ratio = None) -> Generator:
    """Splits a lecture into subreads
    Args:
        seq (str): a DNA sequence
        read_size (int): size of splits
        max_sampling (int): maximum number of samples inside lecture
        shift_ratio (float): ratio of read_size for step in crop
    Raises:
        ValueError: if read is too short
    Yields:
        Generator: subreads collection
    """
    if len(seq) < read_size:
        raise ValueError("Read is too short.")
    if shift_ratio is not None and max_sampling is not None:
        raise ValueError("Provide either shift_ratio or max_sampling, not both.")

    # print("shift_ratio", shift_ratio)
    if shift_ratio is not None:
        assert max_sampling is None
        shift = int(shift_ratio * read_size)
        nb_win = int((len(seq) - read_size) / shift)
    else:
        assert max_sampling is not None
        shift = int((len(seq)-read_size)/max_sampling)
        nb_win = max_sampling

    if shift <= 0 or nb_win<=0:
        raise ValueError(f"Calculated shift {shift} must be positive.")

    # print("shift", shift)
    for i in range(nb_win):
        yield seq[shift*i:shift*i+read_size]


def pattern_filter(substring: str, pattern: list[int]) -> str:
    """Applies a positional filter over a string
    Args:
        substring (str): substring to clean
        pattern (list): integers to be multiplied by
    Returns:
        str: a cleaned kmer
    """
    return ''.join([char * pattern[i] for i, char in enumerate(substring)])


def counter_kmer(entry: str, kmer_size: int, pattern: list[int]) -> Counter:
    """Counts all kmers and filter non-needed ones
    Args:
        entry (str): a subread
        kmer_size (int): k size
        pattern (list[int]): 110110... pattern, to select specific chars in kmer
    Returns:
        Counter: counts of kmers inside subread
    """
    # Defining custom complementarity
    complements: dict = {'A': 'T',
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

    all_kmers: Generator = (entry[i:i+len(pattern)] for i in range(len(entry)-len(pattern)-1))
    counts: Counter = Counter(all_kmers)
    rev_counts: Counter = Counter({revcomp(k, compl=complements): v for k, v in counts.items()})
    counts += rev_counts
    del rev_counts
    if not all(pattern):
        # All positions in pattern should not be kept, we apply filter
        counts = Counter({pattern_filter(k, pattern): v for k, v in counts.items()})
    for filtered_kmer in (alpha * kmer_size for alpha in ['A', 'T', 'C', 'G']):
        if filtered_kmer in counts:
            del counts[filtered_kmer]
    # We treat cases where sequence alphabet is not ATCG
    counts_purged = {}

    for key, count in counts.items():
        list_of_keys = list()
        splitted_key: list = [*key]
        for x in splitted_key:
            if x in ['A', 'T', 'C', 'G']:
                nuct: list = [x]
            elif x == 'U':
                nuct: list = ['T']
            elif x == 'R':
                nuct: list = ['G', 'A']
            elif x == 'Y':
                nuct: list = ['C', 'T']
            elif x == 'K':
                nuct: list = ['G', 'T']
            elif x == 'M':
                nuct: list = ['A', 'C']
            elif x == 'S':
                nuct: list = ['G', 'C']
            elif x == 'W':
                nuct: list = ['A', 'T']
            elif x == 'B':
                nuct: list = ['G', 'T', 'C']
            elif x == 'D':
                nuct: list = ['G', 'T', 'A']
            elif x == 'H':
                nuct: list = ['A', 'T', 'C']
            elif x == 'V':
                nuct: list = ['G', 'A', 'C']
            else:
                nuct: list = ['A', 'T', 'C', 'G']
            # Adding to the keys
            # If list empty
            if len(list_of_keys) == 0:
                list_of_keys = nuct
            else:
                list_of_keys = [new_key+n for n in nuct for new_key in list_of_keys]
        # Updating counts to stay with only ATGC counts
        for prob_key in list_of_keys:
            kmer_number: int = count//len(list_of_keys)
            if prob_key in counts_purged:
                counts_purged[prob_key] += kmer_number
            else:
                counts_purged[prob_key] = kmer_number
    del counts
    # We divide count by the number of keys we end up with to normalize
    return Counter(counts_purged)


@dataclass
class Taxonomy:
    "Modelizes a taxa level"
    code: int | None
    level: str
    name: str
    model_path: str | None
    config_path: str | None

def check_parameters(params: dict) :
    """Lists all conditions where a set of parameters is valid, and accepts the creation if so"""
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
    """Generates a dict of codes for kmers
    Args:
        ksize (int): length of kmer
    Returns:
        dict: kmer:code
    """
    res = {code: encode_kmer(code) for code in map(''.join, product('ATCG', repeat=ksize))}
    return res


def encode_kmer(kmer: str) -> int:
    """Encodes a kmer into base 4 format
    Args:
        kmer (str): a k-sized word composed of A,T,C,G
    Returns:
        int: Encoding of kmer
    """
    mapper: dict = {'A': "0",'C': "1", 'G': "2", 'T': "3",}
    return int(''.join([mapper[k] for k in kmer]))


def taxonomy_information(genome_path: str, tree_struct: Tree) -> tuple[dict, Tree]:
    """Returns taxonomy position information"""
    taxa_family = ['root'] + TAXO_LEVELS
    only_name_file = os.path.splitext(os.path.basename(genome_path))[0]
    # example 'Bacteria_Campylobacterota_Epsilonproteobacteria_Campylobacterales_Campylobacter_hyointestinalis'
    taxo_info: list = only_name_file.split('_')[:5] # only 5 first
    # ['Bacteria', 'Campylobacterota', 'Epsilonproteobacteria', 'Campylobacterales', 'Campylobacter']
    parents = ['Root'] + taxo_info

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

    return dict(zip(TAXO_LEVELS, taxo_info)), tree_struct



def revcomp(string: str, compl=None) -> str:
    """Tries to compute the reverse complement of a sequence
    Args:
        string (str): original character set
        compl (dict, optional): dict of correspondences. Defaults to {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.
    Raises:
        IndexError: Happens if revcomp encounters a char that is not in the dict
    Returns:
        str: the reverse-complemented string
    """
    if compl is None:
        compl = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    try:
        result = ''.join([compl[s] for s in string][::-1])
    except IndexError as exc:
        raise IndexError("Complementarity does not include all chars in sequence.") from exc
    return result
