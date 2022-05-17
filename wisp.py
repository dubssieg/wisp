"This scripts aims to call for command line WISP over a set of files"

###########################################################################################
# This is your main interface. Or rather the call path to your main interface.
# If you want to change any parameter as of database name and such
# head over constants.py : all program constants are stored here.
# Do not change something you're not sure about! If you did so, go to github
# to get a clean copy of the constants file.
# Model parameters and such are stored inside a .json file, "wisp_params" by default
# If you want to change those, please refer to the README.md file!
###########################################################################################

from constants import DATABASE, PARAMS, SAMPLE_PATH, PREFIX_JOB
from os import listdir, system, rename
from python_tools import my_output_msg, my_function_timer, my_logs_clear

# constants ; change those to select database and such


@my_function_timer("Running WISP on unk folder")
def core_call():
    """
    Calls the building and prediction functions with global constants defined above
    If a job fails, skips to the next one
    """
    file_list: list[str] = listdir(
        SAMPLE_PATH)
    # we clean log before entering loop, might be enormous
    my_logs_clear("LOG_wisp.log")

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
