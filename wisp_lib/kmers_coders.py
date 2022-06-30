"Functions to encode and decode kmers"

from functools import cache
from statistics import stdev
from typing import Generator
from collections import Counter
from itertools import product


def recode_kmer_4(input_kmer: str, target_len: int) -> str:
    """Recodes a kmer form its int code, processing at desired size

    Args:
        input_kmer (str): kmer to encode
        target_len (int): len to return

    Returns:
        str: _description_
    """
    if len(input_kmer) == target_len:
        return decode_kmer_4(input_kmer)
    else:
        temp = decode_kmer_4(input_kmer)
        while len(temp) < target_len:
            temp = f"A{temp}"
        return temp


@cache
def encode_kmer_4(kmer: str) -> int:
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


@cache
def decode_kmer_4(kmer: str) -> str:
    """Decodes a kmer from our internal base4 format

    Args:
        kmer (str): a kmer code

    Returns:
        str: the A,T,C,G representation of the kmer
    """
    mapper: dict = {
        '0': "A",
        '3': "T",
        '1': "C",
        '2': "G"
    }
    return ''.join([mapper[k] for k in kmer])


def apply_filter(substring: str, pattern: str) -> str:
    """Set a filter on kmer

    Args:
        substring (str): the kmer to apply filter on
        pattern (str): pattern to multiply by

    Raises:
        ValueError: if kmer is too small for pattern

    Returns:
        str: filtered kmer
    """
    if len(pattern) > len(substring):
        raise ValueError("Substring is too small to apply filter.")
    elif len(substring) == pattern.count('1'):
        return substring
    else:
        return ''.join([substring[i] for i in range(len(substring)) if pattern[i] == '1'])


def reverse_comp(seq: str) -> str:
    """Given a sequence, computes its reversecomp

    Args:
        seq (str): origin sequence

    Returns:
        str: reversed complemented sequence
    """
    try:
        cpl: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([cpl[base]for base in reversed(seq)])
    except Exception as exc:
        raise exc


def read_and_its_compl(entry: str, kmer_size: int, pattern: str) -> Counter:
    """From a read, computes its reverse complement and counts both kmers on read and its reverse

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


def kmer_indexing_brut(entry: str, kmer_size: int) -> Counter:
    """Bruteforce kmer encoding. Used for plots.

    Args:
        entry (str): a read
        kmer_size (int): kmer length

    Returns:
        Counter: kmer counts
    """
    return Counter([entry[k:k+kmer_size] for k in range(len(entry) - kmer_size - 1)])


def encoder_list():
    # maybe this can help gain speed ? specific to k=4
    return {f"{a}{b}{c}{d}": [] for a in ['A', 'T', 'G', 'C'] for b in ['A', 'T', 'G', 'C'] for c in ['A', 'T', 'G', 'C']for d in ['A', 'T', 'G', 'C']}


def kmer_indexing_10000(entry: str, kmer_size: int) -> dict:
    """Bruteforce kmer encoding. Used for plots.

    Args:
        entry (str): a read
        kmer_size (int): kmer length

    Returns:
        Counter: kmer counts
    """
    entries = optimal_splitting(entry, 10000, 100)
    entrs = encoder_list()
    for entr in entries:
        counts = Counter([entr[k:k+kmer_size]
                         for k in range(len(entr) - kmer_size - 1)])
        for k, v in counts.items():
            entrs[k].append(float(v))
    return {key: stdev(value) for key, value in entrs.items()}


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


def counter_fast(entry: str, kmer_size: int, pattern: str) -> Counter:
    """Computes kmer counts in read

    Args:
        entry (str): a read
        kmer_size (int): size of words we use
        pattern (str): pattern to apply to kmers

    Returns:
        Counter: counts of kmers
    """
    all_kmers: list = [entry[i:i+len(pattern)]
                       for i in range(len(entry)-len(pattern)-1)]
    if pattern != len(pattern) * '1':
        all_kmers = [apply_filter(kmer, pattern) for kmer in all_kmers]
    filtering = [alpha * kmer_size for alpha in ['A', 'T', 'C', 'G']]
    return Counter(x for x in all_kmers if x not in filtering)


def apply_filter_fast(substring: str, pattern: str) -> str:
    """From a set of chars, filter those which are discriminated by pattern

    Args:
        substring (str): set of chars we need to expurge chars from
        pattern (str): pattern to select chars

    Returns:
        str: epurged str
    """
    # function that takes most of the time
    return ''.join([c for i, c in enumerate(substring) if pattern[i] == '1'])

####################################################################################


def counter_ultrafast(entry: str, kmer_size: int, pattern: list) -> Counter:
    """Counts all kmers and filter non-needed ones

    Args:
        entry (str): a subread
        kmer_size (int): k size
        pattern (str): 110110... pattern, to select specific chars in kmer

    Returns:
        Counter: counts of kmers inside subread
    """
    all_kmers = (entry[i:i+len(pattern)]
                 for i in range(len(entry)-len(pattern)-1))
    counts = Counter(all_kmers)
    counts_reverse = Counter({reverse_comp(k): v for k, v in counts.items()})
    all_counts = counts + counts_reverse
    if pattern != len(pattern) * '1':
        counts = Counter({ultrafast_filter(k, pattern): v for k, v in all_counts.items()})
    for filtered_kmer in (alpha * kmer_size for alpha in ['A', 'T', 'C', 'G']):
        del counts[filtered_kmer]
    return counts


def ultrafast_filter(substring: str, pattern: list) -> str:
    """Applies a positional filter over a string (! NEEDS REFINING)

    Args:
        substring (str): substring to clean
        pattern (list): integers to be multiplied by

    Returns:
        str: a cleaned kmer
    """
    return ''.join([char * pattern[i] for i, char in enumerate(substring)])


def splitting_generator(seq: str, window_size: int, max_sampling: int) -> Generator:
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


def encoder(ksize: int) -> dict:
    """Generates a dict of codes for kmers

    Args:
        ksize (int): length of kmer

    Returns:
        dict: kmer:code
    """
    return {code: encode_kmer_4(code) for code in map(''.join, product('ATCG', repeat=ksize))}
