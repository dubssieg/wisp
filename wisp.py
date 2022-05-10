"Calls for command line wisp over a set of files"

from os import listdir, system, rename
from python_tools import my_output_msg, my_function_timer

DATABASE = "4k_5for10_10000pb"
PARAMS = "wisp_params.json"
JOB = "fasta_for_read"
PATH = "/udd/sidubois/Stage/Genomes/unk/"


def rename_genomes(path):
    "If preprocessing of names is necessary for you"
    files = [file for file in listdir(path)]

    for file in files:
        new_name = file.replace('_', '-')
        rename(f"{path}{file}",
               f"{path}{new_name}")
    return listdir(path)


@my_function_timer("Running WISP on unk folder")
def core_call():
    file_list = listdir(PATH)  # rename_genomes(path) as alternate call

    for i, file in enumerate(file_list):
        system(
            f"python main.py {DATABASE} {PARAMS} 4k_5for10_{file[:-4]} -f {file}")
    my_output_msg(
        f"Sucessfully processed {len(file_list)} genomes. Results are in output/ folder")


if __name__ == "__main__":
    "Executes main procedure"
    core_call()
