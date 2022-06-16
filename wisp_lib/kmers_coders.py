"Functions to encode and decode kmers"

from functools import cache
from typing import Generator
from khmer import Countgraph
from collections import Counter
from time import monotonic
from random import choice
from itertools import product


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


def kmer_indexing_brut(entry: str, kmer_size: int):
    return Counter([entry[k:k+kmer_size] for k in range(len(entry) - kmer_size - 1)])


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


def counter_fast(entry: str, kmer_size: int, pattern: str):
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
        counts = Counter({ultrafast_filter(k, pattern)                         : v for k, v in all_counts.items()})
    for f in (alpha * kmer_size for alpha in ['A', 'T', 'C', 'G']):
        del counts[f]
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

####################################################################################


if __name__ == "__main__":
    seq = ''.join([choice(['A', 'T', 'C', 'G']) for _ in range(1000000)])
    my_encoder = encoder(5)

    base = monotonic()

    splits = splitting_generator(seq, 10000, 500)
    counters = [kmer_indexing_canonical(split, 5, "111011")
                for split in splits]

    print(f"kmer_canonique:generator en : {monotonic()-base}")
    base = monotonic()

    splits = splitting_generator(seq, 10000, 500)
    counters = [counter_ultrafast(split, 5, [1, 1, 1, 0, 1, 1])
                for split in splits]

    print(f"counter_ultrafast:generator en : {monotonic()-base}")
    base = monotonic()

    encoded = [{my_encoder[k]:v for k, v in cts.items()} for cts in counters]

    #print(f"encode:generator en : {monotonic()-base}")
    base = monotonic()

    encoded = [{encode_kmer_4(k): v for k, v in cts.items()}
               for cts in counters]

    #print(f"encode:regular en : {monotonic()-base}")
    base = monotonic()

    splits = optimal_splitting(seq, 10000, 500)
    [counter_ultrafast(split, 5, [1, 1, 1, 0, 1, 1]) for split in splits]

    print(f"counter_ultrafast:regular en : {monotonic()-base}")

    base = monotonic()

    splits = optimal_splitting(seq, 10000, 500)
    [counter_fast(split, 5, '111011') for split in splits]

    print(f"counter_fast en : {monotonic()-base}")

    base = monotonic()

    splits = optimal_splitting(seq, 10000, 500)
    [kmer_indexing_brut(split, 5, '111011') for split in splits]

    print(f"kmer_indexing_brut en : {monotonic()-base}")
