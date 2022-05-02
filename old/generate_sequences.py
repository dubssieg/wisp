import random
from python_tools import my_output_msg

# consts
SIZE_KMERS: int = 100
FAM_NUMBER: int = 4
SIZE_BP: int = 1500
VARIABILITY: float = 0.01
ALT_NUMBER: int = 5000

######################### GENERATE SEQUENCES ###########################


def my_seq_generator(a: int = FAM_NUMBER, x: int = SIZE_BP):
    # generate a base sequence of X bp from which we will do SNP/indels
    return {f"seq{i}": ''.join([random.choice(["A", "T", "C", "G"]) for _ in range(x)]) for i in range(a)}


def alteration(seq_dict: dict = my_seq_generator(), var: float = VARIABILITY, n: int = ALT_NUMBER):
    # generates n sequences by seq in seq_dict with variation of var
    return {f"{key}_{y}": ''.join([char if random.random() > var else random.choice(['', '', '', '', f"{char}A", f"{char}T", f"{char}C", f"{char}G", "A", "T", "C", "G"]) for char in value]) for y in range(n) for key, value in seq_dict.items()}

######################### GENERATE DATASETS ###########################


def generate_dataset(data: dict, job: str, list_sp_codes: dict, ext: str) -> None:
    "data must be dict of dicts"
    my_output_msg("Starting genome mapping...")
    lines: list = []
    for key, value in data.items():
        attributes = ' '.join(
            [f"{subkey}:{subvalue}" for subkey, subvalue in value.items()])
        lines.append(f"{list_sp_codes[key.split('_')[0]]} {attributes}")
    with open(f"data/{job}.txt.{ext}", "w") as writer:
        writer.write('\n'.join(lines))


def save_seqs(seq_dict: dict = alteration()):
    lines = []
    trad: dict = {
        'A': 0,
        'T': 1,
        'C': 2,
        'G': 3
    }
    for key, value in seq_dict.items():
        dic = {f"read{i}": value[100*i:100*i + SIZE_KMERS]
               for i in range(0, 10, 1)}
        val = ' '.join(
            [f"{i}:{trad[dic['read0'][i]]}" for i in range(len(dic['read0']))])
        #val2 = f"0:{dic['read0'][0]}"
        #val3 = f"1:{dic['read1']} 2:{dic['read2']} 3:{dic['read3']} 4:{dic['read4']} 5:{dic['read5']} 6:{dic['read6']} 7:{dic['read7']} 8:{dic['read8']} 9:{dic['read9']}"
        line = f"{key[3]} {val}"  # LIBSVM formating
        lines.append(line)
    random.shuffle(lines)
    test = lines[:int(len(lines)/2)]
    train = lines[int(len(lines)/2):]
    with open("data/generated_seqs.txt.test", "w") as writer:
        writer.write('\n'.join(test))
    with open("data/generated_seqs.txt.train", "w") as writer:
        writer.write('\n'.join(train))


if __name__ == "__main__":
    "Executes main procedure"
    save_seqs()
    print("Job done")
