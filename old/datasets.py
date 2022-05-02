
"""
@my_function_timer("Building datasets")
def make_datasets_with_multithread(input_dir: str, path: str, job_name: str, func, ratio: float, kmer_size: int):
    "Creating datasets with new methods"
    sp_map = mapping_sp(f"{input_dir}train/", path, job_name)
    KMER_SIZE = kmer_size
    for type in ['test', 'train']:
        my_sp = list(call_loader(f"{input_dir}{type}/", kmer_size))
        lines = []
        for sp in my_sp:
            if func != None:
                sp.update_counts(func, ratio)
            lines.append(sp.encoding_mapping(
                sp_map[sp.specie.split('_')[0]]))
        output_writing(lines, path, job_name, type)
        print(
            f"{type} : loaded {len(my_sp)} genomes under {len(sp_map)} distinct clades")
"""

"""
@my_function_timer("Building datasets")
def make_datasets(input_dir: str, path: str, job_name: str, func, ratio, kmer_size):
    "Creating datasets with new methods"
    sp_map = mapping_sp(f"{input_dir}train/", path, job_name)

    for type in ['train', 'test']:
        my_sp = load_list_sample(f"{input_dir}{type}/", kmer_size)
        lines = []
        for sp in my_sp:
            if func != None:
                sp.update_counts(func, ratio)
            lines.append(sp.encoding_mapping(sp_map[sp.specie.split('_')[0]]))
        output_writing(lines, path, job_name, type)
        print(
            f"{type} : loaded {len(my_sp)} genomes under {len(sp_map)} distinct clades")
"""

"""
@my_function_timer("Building datasets")
def make_datasets_with_sampling(input_dir: str, path: str, job_name: str, func, ratio, kmer_size):
    "Creating datasets with new methods"
    sp_map = mapping_sp(f"{input_dir}train/", path, job_name)

    for type in ['train']:
        my_sp = list(call_loader(f"{input_dir}{type}/", kmer_size))
        lines = []
        for sp in my_sp:
            print(sp)
            if func != None:
                sp.update_counts(func, ratio)
            lines.append(sp.encoding_mapping(sp_map[sp.specie.split('_')[0]]))
        output_writing(lines, path, job_name, type)
        print(
            f"{type} : loaded {len(my_sp)} genomes under {len(sp_map)} distinct clades")

    for type in ['test']:
        my_sp = load_list_sample(f"{input_dir}{type}/", kmer_size)
        lines = []
        for sp in my_sp:
            my_ssp = generate_diversity(sp, 20, 5, 10000)
            if func != None:
                sp.update_counts(func, ratio)
            for s in my_ssp:
                lines.append(s.encoding_mapping(
                    sp_map[s.specie.split('_')[0]]))
        output_writing(lines, path, job_name, type)
        print(
            f"{type} : loaded {len(my_sp)} genomes under {len(sp_map)} distinct clades")
"""
