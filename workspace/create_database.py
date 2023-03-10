"Creates a json database"
from collections import Counter
from json import load, dump
from os import path
from pathlib import Path
from typing import Generator
from itertools import product
from Bio import SeqIO
from tharospytools import revcomp


def validate_parameters(params: dict) -> bool:
    "Lists all conditions where a set of parameters is valid, and accepts the creation if so"
    return all(
        [
            # verifies that the pattern length respects ksize
            sum(params['pattern']) == params['ksize'],
        ]
    )


def encoder(ksize: int) -> dict:
    """Generates a dict of codes for kmers

    Args:
        ksize (int): length of kmer

    Returns:
        dict: kmer:code
    """
    return {code: encode_kmer(code) for code in map(''.join, product('ATCG', repeat=ksize))}


def encode_kmer(kmer: str) -> int:
    """Encodes a kmer into base 4 format

    Args:
        kmer (str): a k-sized word composed of A,T,C,G

    Returns:
        int: Encoding of kmer
    """
    mapper: dict = {
        'A': "0",
        'T': "3",
        'C': "1",
        'G': "2"
    }
    return int(''.join([mapper[k] for k in kmer]))


def taxonomy_information(genome_path: str) -> dict:
    "Returns taxonomy position information"
    taxo_info: list = Path(genome_path).stem.split('_')
    return {
        'domain': taxo_info[0],
        'phylum': taxo_info[1],
        'group': taxo_info[2],
        'order': taxo_info[3],
        'family': taxo_info[4]
    }


def build_database(params_file: str, database_name: str, input_data: list[str]) -> str:
    "Builds a json file with taxa levels as dict information"
    # Loading params file
    with open(params_file, 'r', encoding='utf-8') as pfile:
        params: dict = load(pfile)

    # Guard to check if params are acceptable
    if not validate_parameters(params):
        raise RuntimeError("Incorrect parameter file")

    # creating encoder
    my_encoder: dict = encoder(ksize=params['ksize'])

    # Writing the database
    json_datas: list = list()
    Path(f"{path.dirname(__file__)}/databases/").mkdir(parents=True, exist_ok=True)
    with open(output_path := f'{path.dirname(__file__)}/databases/{database_name}.json', 'w', encoding='utf-8') as jdb:

        # iterating over input genomes
        for genome in input_data:
            with open(genome, 'r', encoding='utf-8') as freader:
                genome_data: dict = {fasta.id: str(fasta.seq)
                                     for fasta in SeqIO.parse(freader, 'fasta')}
            for id_sequence, dna_sequence in genome_data.items():
                # Splitting of reads
                if len(dna_sequence) >= params['read_size']:
                    all_reads = splitting(
                        dna_sequence,
                        params['read_size'],
                        params['sampling']
                    )
                    # Counting kmers inside each read
                    counters: list[Counter] = [
                        counter(
                            read,
                            params['ksize'],
                            params['pattern']
                        ) for read in all_reads
                    ]
                    del all_reads

                    # Encoding reads for XGBoost
                    encoded: list = [{my_encoder[k]:v for k, v in cts.items()}
                                     for cts in counters]
                    del counters

                    # Dumping in output file
                    json_datas.append(
                        {**taxonomy_information(genome), 'datas': encoded})
                    del encoded
        dump({"mappings": mapping_sp(json_datas), "datas": json_datas}, jdb)
    return output_path


def mapping_sp(datas: list[dict]) -> dict:
    """Given a dataset, creates numbers for each level of classification

    Args:
        datas (list[dict]): a set of assigned reads

    Returns:
        dict: a grouped-by-level list of codes
    """
    taxa: list[str] = [
        'domain',
        'phylum',
        'group',
        'order',
        'family'
    ]
    taxa_codes: dict = {taxon: {'number_taxa': 0} for taxon in taxa}
    for sample in datas:
        for key, value in sample.items():
            if key in taxa:
                if not value in (taxa_level := taxa_codes[key]):
                    taxa_level[value] = taxa_level['number_taxa']
                    taxa_level['number_taxa'] += 1
    return taxa_codes


def splitting(seq: str, window_size: int, max_sampling: int) -> Generator:
    """Splits a lecture into subreads

    Args:
        seq (str): a DNA sequence
        window_size (int): size of splits
        max_sampling (int): maximum number of samples inside lecture

    Raises:
        ValueError: if read is too short

    Yields:
        Generator: subreads collection
    """
    if len(seq) < window_size:
        raise ValueError("Read is too short.")
    shift: int = int((len(seq)-window_size)/max_sampling)
    for i in range(max_sampling):
        yield seq[shift*i:shift*i+window_size]


def pattern_filter(substring: str, pattern: list[int]) -> str:
    """Applies a positional filter over a string

    Args:
        substring (str): substring to clean
        pattern (list): integers to be multiplied by

    Returns:
        str: a cleaned kmer
    """
    return ''.join([char * pattern[i] for i, char in enumerate(substring)])


def counter(entry: str, kmer_size: int, pattern: list[int]) -> Counter:
    """Counts all kmers and filter non-needed ones

    Args:
        entry (str): a subread
        kmer_size (int): k size
        pattern (list[int]): 110110... pattern, to select specific chars in kmer

    Returns:
        Counter: counts of kmers inside subread
    """

    all_kmers: Generator = (entry[i:i+len(pattern)]
                            for i in range(len(entry)-len(pattern)-1))
    counts: Counter = Counter(all_kmers)
    rev_counts: Counter = Counter({revcomp(k): v for k, v in counts.items()})
    counts += rev_counts
    del rev_counts
    if not all(pattern):
        # All positions in pattern should not be kept, we apply filter
        counts = Counter({pattern_filter(k, pattern): v for k, v in counts.items()})
    for filtered_kmer in (alpha * kmer_size for alpha in ['A', 'T', 'C', 'G']):
        del counts[filtered_kmer]
    return counts
