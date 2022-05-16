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


def kmer_indexing(entry: str, kmer_size: int):
    ksize = kmer_size
    nkmers = 4**ksize
    tablesize = nkmers + 10
    cg = Countgraph(ksize, tablesize, 1)
    for k in range(len(entry) - ksize + 1):
        cg.count(entry[k:k+ksize])
    return Counter({cg.reverse_hash(i): cg.get(i) for i in range(nkmers) if cg.get(i)})
