#!/usr/bin/env python3.10
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

from constants import DATABASE, PARAMS, SAMPLE_PATH, PREFIX_JOB, LEVELS
from os import listdir
from sys import executable
from python_tools import my_output_msg, my_function_timer, my_logs_global_config, my_futures_collector
from argparse import ArgumentParser
from datetime import datetime
import shlex
import subprocess


# constants ; change those to select database and such


@my_function_timer("WISP run")
def core_call(multithreading_state: int, building_state: bool) -> None:
    """
    Calls the building and prediction functions with global constants defined above
    If a job fails, skips to the next one
    """
    try:
        match building_state:
            case True:
                communicators = my_futures_collector(subprocess.Popen, [
                    [shlex.split(f"{executable} force_build.py {DATABASE} {PARAMS} -l {[level]}")] for level in LEVELS], multithreading_state)
            case False:
                file_list: list[str] = listdir(SAMPLE_PATH)
                communicators = my_futures_collector(subprocess.Popen, [[
                    shlex.split(f"{executable} main.py {DATABASE} {PARAMS} {PREFIX_JOB}_{file[:-4]} -f {file}")] for file in file_list], multithreading_state)
        retcodes = [p.communicate() for p in communicators]
    except Exception as exc:
        raise BaseException("Job failed") from exc
    finally:
        for i, code in enumerate(retcodes):
            my_output_msg(f"Job {i} : {code}")


if __name__ == "__main__":
    "Executes main procedure"
    # we clean log before entering loop, might be enormous
    parser = ArgumentParser()
    parser.add_argument(
        "-v", "--verbose", help="Detailed log file", action='store_true')
    parser.add_argument(
        "-b", "--build", help="Calls for database building instead of prediciton", action='store_true')
    parser.add_argument(
        "-t", "--multithreading", type=int, default=1, help="Gives a thread number for WISP")
    args = parser.parse_args()
    my_logs_global_config("LOG_wisp", verbose=args.verbose)
    core_call(args.multithreading, args.build)
