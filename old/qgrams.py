import collections
import generate_sequences
import json
from statistics import mean
import random
from python_tools import my_output_msg
import functools

KMER_SIZE: int = 4


@functools.cache
def encode_kmer_4(kmer: str) -> int:
    "Encodes a kmer into base 4 format"
    mapper: dict = {
        'A': "0",
        'T': "3",
        'C': "1",
        'G': "2"
    }
    return int(''.join([mapper[k] for k in kmer]))


@functools.cache
def decode_kmer_4(kmer: str) -> str:
    "Decodes a kmer from base 4 format"
    mapper: dict = {
        '0': "A",
        '3': "T",
        '1': "C",
        '2': "G"
    }
    return ''.join([mapper[k] for k in kmer])


def indexing(df: dict):
    """
    Returns a dict of dicts containing for each species its encoded kmars and their numbers

    * df must be a dict of sequences
    """
    all_reads: dict = {}
    for key, value in df.items():
        read = f"{value}{value[-(KMER_SIZE-1):]}"
        encode_kmers = [encode_kmer(read[i:i+KMER_SIZE])
                        for i in range(len(value))]
        sum_kmers = collections.Counter(encode_kmers)
        all_reads[key] = {f"{k}": f"{v}" for k, v in sum_kmers.items() if int(
            v) > mean(sum_kmers.values())*(1.5)}
    return all_reads


def indexing_by_signature(df: dict, func=mean, ratio: float = 1.5):
    """
    Returns a dict of dicts containing for each species its encoded kmars and their numbers

    * df must be a dict of sequences
    """
    all_reads: dict = {}
    for key, value in df.items():
        read: str = f"{value}{value[-(KMER_SIZE-1):]}"
        encode_kmers: list = [encode_kmer(
            read[i:i+KMER_SIZE]) for i in range(len(value))]
        sum_kmers = collections.Counter(encode_kmers)
        all_reads[key] = {f"{k}": f"{round(float((v)/len(read)*100000),2)}" for k,
                          v in sum_kmers.items() if int(v) > func(sum_kmers.values())*(ratio)}
    return all_reads


def indexing_by_signature_with_subsampling(df: dict, func=mean, ratio: float = 1.5, sample_ratio: float = 0.0005, size_reads: int = 10000):
    """
    Returns a dict of dicts containing for each species its encoded kmars and their numbers

    * df must be a dict of sequences
    * func must be an agregation func
    * ratio is float and is a multiplier for func
    * sample_ratio gives number of reads subsampled out from genome size
    * size_reads is the size in bp for each individual read
    """
    all_reads: dict = {}
    for key, value in df.items():
        read: str = f"{value}{value[-(KMER_SIZE-1):]}"
        number_samples: int = int(len(read)*sample_ratio)
        for sample in range(number_samples):
            x = random.randrange(0, len(read)-size_reads-KMER_SIZE)
            subread = read[x:x+size_reads+KMER_SIZE]
            encode_kmers: list = [encode_kmer(
                subread[i:i+KMER_SIZE]) for i in range(len(subread))]
            sum_kmers = collections.Counter(encode_kmers)
            all_reads[f"{key}_{sample}"] = {
                f"{k}": f"{round(float((v)/len(subread)*100000),2)}" for k, v in sum_kmers.items() if int(v) > func(sum_kmers.values())*(ratio)}
    return all_reads


def indexing_by_signature_with_subsampling_linear(df: dict, func=mean, ratio: float = 1.3, sample_ratio: float = 0.00005, size_reads: int = 5000):
    """
    Returns a dict of dicts containing for each species its encoded kmars and their numbers

    * df must be a dict of sequences
    * func must be an agregation func
    * ratio is float and is a multiplier for func
    * sample_ratio gives number of reads subsampled out from genome size
    * size_reads is the size in bp for each individual read
    """
    all_reads: dict = {}
    final_read: dict = {}
    for key, value in df.items():
        read: str = f"{value}{value[-(KMER_SIZE-1):]}"
        number_samples: int = int(len(read)*sample_ratio)
        for sample in range(number_samples):
            x = random.randrange(0, len(read)-size_reads-KMER_SIZE)
            subread = read[x:x+size_reads+KMER_SIZE]
            encode_kmers: list = [encode_kmer_4(
                subread[i:i+KMER_SIZE]) for i in range(len(subread))]
            sum_kmers = collections.Counter(encode_kmers)
            all_reads[f"{key}_{sample}"] = {f"{k}": f"{all_reads.get(key,{k:0}).get(k,0) + round(float((v)/len(subread)*100000),2)}" for k, v in sum_kmers.items(
            ) if int(v) > func(sum_kmers.values())*(1.5) or int(v) < func(sum_kmers.values())*(0.5)}
        #final_read[key] = {k: v for ld in list(all_reads.items()) for k, v in ld.items()}
        ret = {}
        for _, value in all_reads.items():
            ret = ret | value
        final_read[key] = ret
        #dict(functools.reduce(operator.add, map(collections.Counter, all_reads)))
    return final_read


def main(index_method=indexing_by_signature_with_subsampling, with_test: bool = False):
    my_output_msg("Rebuilding datasets and indexes...")

    # genra mapping
    list_sp = load_sequences.species_map("/udd/sidubois/Stage/Genomes/")
    file_manager = open("/udd/sidubois/Stage/wisp/data/mapping.json", "w")
    file_manager = json.dump(list_sp, file_manager)

    # training set
    df = load_sequences.sequences("/udd/sidubois/Stage/Genomes/")
    generate_sequences.generate_dataset(
        index_method(df), "two_bacteries", list_sp, 'train')

    if(with_test):
        # validation set
        df_generated = load_sequences.sequences(
            "/udd/sidubois/Stage/Genomes_test/")
        generate_sequences.generate_dataset(index_method(
            df_generated), "two_bacteries", list_sp, 'test')

    my_output_msg("Job done")


if __name__ == "__main__":
    main()
