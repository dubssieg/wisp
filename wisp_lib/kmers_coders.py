"Functions to encode and decode kmers"

from functools import cache
from khmer import Countgraph
from collections import Counter


def recode_kmer_4(input: str, target_len: int):
    if len(input) == target_len:
        return decode_kmer_4(input)
    else:
        temp = decode_kmer_4(input)
        while len(temp) < target_len:
            temp = f"A{temp}"
        return temp


@cache
def encode_kmer_4(kmer: str) -> int:
    "Encodes a kmer into base 4 format"
    mapper: dict = {
        'A': "0",
        'T': "3",
        'C': "1",
        'G': "2"
    }
    return int(''.join([mapper[k] for k in kmer]))


@cache
def decode_kmer_4(kmer: str) -> str:
    "Decodes a kmer from base 4 format"
    mapper: dict = {
        '0': "A",
        '3': "T",
        '1': "C",
        '2': "G"
    }
    return ''.join([mapper[k] for k in kmer])


def my_encoder_k4():
    # maybe this can help gain speed ? specific to k=4
    return {f"{a}{b}{c}{d}": encode_kmer_4(f"{a}{b}{c}{d}") for a in ['A', 'T', 'G', 'C'] for b in ['A', 'T', 'G', 'C'] for c in ['A', 'T', 'G', 'C']for d in ['A', 'T', 'G', 'C']}


def apply_filter(substring: str, pattern: str) -> str:
    if len(pattern) > len(substring):
        raise ValueError("Substring is too small to apply filter.")
    elif len(substring) == pattern.count('1'):
        return substring
    else:
        return ''.join([substring[i] for i in range(len(substring)) if pattern[i] == '1'])


def reverse_comp(seq: str) -> str:
    cpl: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([cpl[base]for base in reversed(seq)])  # reverse complement


def read_and_its_compl(entry: str, kmer_size: int, pattern: str) -> Counter:
    """Froma read, computes its reverse complement and counts both kmers on read and its reverse

    Args:
        entry (str): a read to be computed
        kmer_size (int): size of window we're reading with
        pattern (str): a seed in format 1:keep and 0:ignore

    Raises:
        ValueError: If filter is not set accordingly to ksize, we raise an error

    Returns:
        Counter: sum of kmers from read and its reverse comp
    """
    if pattern.count('1') != kmer_size:
        raise ValueError("Filter does not match ksize.")
    else:
        entry_reverse = reverse_comp(entry)
        return Counter([apply_filter(entry[k:k+len(pattern)], pattern) for k in range(len(entry) - len(pattern) - 1) if apply_filter(entry[k:k+len(pattern)], pattern).isalpha() and not apply_filter(entry[k:k+len(pattern)], pattern) == len(apply_filter(entry[k:k+len(pattern)], pattern)) * apply_filter(entry[k:k+len(pattern)], pattern)[0]]) + Counter([apply_filter(entry_reverse[k:k+len(pattern)], pattern) for k in range(len(entry_reverse) - len(pattern) - 1) if apply_filter(entry_reverse[k:k+len(pattern)], pattern).isalpha() and not apply_filter(entry_reverse[k:k+len(pattern)], pattern) == len(apply_filter(entry_reverse[k:k+len(pattern)], pattern)) * apply_filter(entry_reverse[k:k+len(pattern)], pattern)[0]])


def kmer_2soluces(entry: str, kmer_size: int, pattern: str, inverted: bool = False):
    """Depending of state of bool, counts kmers in 3'->5' or in 5'->3'

    Args:
        entry (str): a read to be computed
        kmer_size (int): size of window we're reading with
        pattern (str): a seed in format 1:keep and 0:ignore
        inverted (bool, optional): reverse sequence beforehand. Defaults to False.

    Raises:
        ValueError: If filter is not set accordingly to ksize, we raise an error

    Returns:
        Counter: sum of kmers from read or its reverse comp
    """
    if pattern.count('1') != kmer_size:
        raise ValueError("Filter does not match ksize.")
    else:
        if inverted:
            entry = reverse_comp(entry)
        return Counter([apply_filter(entry[k:k+len(pattern)], pattern) for k in range(len(entry) - len(pattern) - 1) if apply_filter(entry[k:k+len(pattern)], pattern).isalpha() and not apply_filter(entry[k:k+len(pattern)], pattern) == len(apply_filter(entry[k:k+len(pattern)], pattern)) * apply_filter(entry[k:k+len(pattern)], pattern)[0]])


def kmer_indexing_canonical(entry: str, kmer_size: int, pattern: str) -> Counter:
    if pattern.count('1') != kmer_size:
        raise ValueError("Filter does not match ksize.")
    else:
        ksize = kmer_size
        nkmers = 4**ksize
        tablesize = nkmers + 10
        cg = Countgraph(ksize, tablesize, 1)
        for k in range(len(entry) - len(pattern) - 1):
            my_kmer = apply_filter(entry[k:k+len(pattern)], pattern)
            if my_kmer.isalpha() and not my_kmer == len(my_kmer) * my_kmer[0]:
                # fixes rare issue where a \n was integrated
                # and checks for mononucleotid patterns
                cg.count(my_kmer)
        return Counter({cg.reverse_hash(i): cg.get(i) for i in range(nkmers) if cg.get(i)})


def kmer_indexing_brut(entry: str, kmer_size: int, pattern: str):
    if pattern.count('1') != kmer_size:
        raise ValueError("Filter does not match ksize.")
    else:
        return Counter([apply_filter(entry[k:k+len(pattern)], pattern) for k in range(len(entry) - len(pattern) - 1) if apply_filter(entry[k:k+len(pattern)], pattern).isalpha() and not apply_filter(entry[k:k+len(pattern)], pattern) == len(apply_filter(entry[k:k+len(pattern)], pattern)) * apply_filter(entry[k:k+len(pattern)], pattern)[0]])


def optimal_splitting(seq: str, window_size: int, max_sampling: int) -> set[str]:
    """_summary_

    Args:
        seq (str): _description_
        window_size (int): _description_
        max_sampling (int): _description_

    Raises:
        ValueError: _description_

    Returns:
        set[str]: unique reads covering sequence
    """
    if len(seq) < window_size:
        raise ValueError("Read is too short.")
    shift: int = int((len(seq)-window_size)/max_sampling)
    return set([seq[shift*i:shift*i+window_size] for i in range(max_sampling)])
