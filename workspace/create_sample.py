"Creates the sample dataset to be predicted afterwards"
from time import time
from itertools import product
from typing import Generator
from pathlib import Path
from os import path
from json import load
from collections import Counter
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


def build_sample(params_file: str, dna_sequence: str, id_sequence: str) -> str:
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
    Path(f"{path.dirname(__file__)}/databases/").mkdir(parents=True, exist_ok=True)
    with open(output_path := f"{path.dirname(__file__)}/databases/unk_sample_{str(time()).replace('.','_')}_{id_sequence.replace(' ','_')}.txt", 'w', encoding='utf-8') as jdb:

        all_reads = splitting(
            dna_sequence.upper(),
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

        for sample in encoded:
            # Each read is a dict with code:count for kmer
            jdb.write(
                f"0 {' '.join([str(k)+':'+str(v) for k,v in sample.items()])} #{id_sequence}\n")

        del encoded

    return output_path


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
    # Defining custom complementarity
    complements: dict = {
        'A': 'T',
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
        'N': 'N'
    }

    all_kmers: Generator = (entry[i:i+len(pattern)]
                            for i in range(len(entry)-len(pattern)-1))
    counts: Counter = Counter(all_kmers)
    rev_counts: Counter = Counter(
        {revcomp(k, compl=complements): v for k, v in counts.items()})
    counts += rev_counts
    del rev_counts
    if not all(pattern):
        # All positions in pattern should not be kept, we apply filter
        counts = Counter({pattern_filter(k, pattern): v for k, v in counts.items()})
    for filtered_kmer in (alpha * kmer_size for alpha in ['A', 'T', 'C', 'G']):
        if filtered_kmer in counts:
            del counts[filtered_kmer]
    # We treat cases where sequence alphabet is not ATCG
    counts_purged: dict = {}
    for key, count in counts.items():
        list_of_keys: list = list()
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
                list_of_keys = [
                    new_key+n for n in nuct for new_key in list_of_keys]
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
