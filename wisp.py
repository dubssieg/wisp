"Calls for command line wisp over a set of files"

from os import listdir, system, rename
from python_tools import my_output_msg, my_function_timer, my_logs_clear

# constants ; change those to select database and such
DATABASE: str = "my_db"
PARAMS: str = "wisp_params.json"
PATH: str = "/udd/sidubois/Stage/Genomes/unk/"
PREFIX_JOB: str = "plot"


def rename_genomes(path: str):
    "If preprocessing of names is necessary for you"
    files: list[str] = [file for file in listdir(path)]

    for file in files:
        new_name = file.replace('_', '-')
        rename(f"{path}{file}",
               f"{path}{new_name}")
    return listdir(path)


@my_function_timer("Running WISP on unk folder")
def core_call():
    """
    Calls the building and prediction functions with global constants defined above
    If a job fails, skips to the next one
    """
    file_list: list[str] = listdir(
        PATH)  # rename_genomes(path) as alternate call
    my_logs_clear("LOG_wisp.log")  # we clean before entering loop

    for i, file in enumerate(file_list):
        try:  # if a job happens to fail, you can check the .log file to check the crash cause
            system(
                f"python main.py {DATABASE} {PARAMS} {PREFIX_JOB}_{file[:-4]} -f {file}")
            my_output_msg(
                f"Sucessfully processed {len(file_list)} genomes. Results are in output/ folder")
        except:
            my_output_msg(f"Job failed for sample number {i} ({file})")


if __name__ == "__main__":
    "Executes main procedure"
    core_call()
