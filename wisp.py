#!/usr/bin/env python3.10
"This scripts aims to call for command line WISP over a set of files"
###########################################################################################
# This is your main interface. Or rather the call path to your main interface.
# Model parameters and such are stored inside a .json file, "wisp_params" by default
# If you want to change those, please refer to the README.md file!
###########################################################################################

from os import listdir
from sys import executable
from python_tools import my_output_msg, my_function_timer, my_logs_global_config, my_futures_collector
from wisp_lib import load_json
from argparse import ArgumentParser
from traceback import format_exc
import shlex
import subprocess


# constants ; change those to select database and such


@my_function_timer("WISP run")
def core_call(exclusion: str, multithreading_state: int, building_state: bool, params: str, taxas_levels: list[str], db: str, unk_path: str, job_prefix: str) -> None:
    """
    Calls the building and prediction functions with global constants defined above
    If a job fails, skips to the next one
    """
    try:
        match building_state:
            case True:
                communicators = my_futures_collector(subprocess.Popen, [
                    [shlex.split(f"{executable} force_build.py {db} {params} -l {[level]}")] for level in taxas_levels], multithreading_state)
            case False:
                file_list: list[str] = listdir(unk_path)
                communicators = my_futures_collector(subprocess.Popen, [[
                    shlex.split(f"{executable} main.py {db} {params} {job_prefix}_{file[:-4]} -f {file} -e {exclusion}")] for file in file_list], multithreading_state)
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
        "params", help="path to a params .json file", type=str)
    parser.add_argument(
        "-v", "--verbose", help="Detailed log file", action='store_true')
    parser.add_argument(
        "-b", "--build", help="Calls for database building instead of prediciton", action='store_true')
    parser.add_argument(
        "-t", "--multithreading", type=int, default=1, help="Gives a thread number for WISP")
    parser.add_argument(
        "-e", "--exclude", type=str, help="Used for one-vs-all tests")
    args = parser.parse_args()

    try:
        my_params: dict = load_json(args.params)

        LOG_FILE: str = my_params['log_file']
        DATABASE: str = my_params['db_name']
        TAXAS_LEVELS: list[str] = my_params['levels_list']
        SAMPLES_PATH: str = my_params['input_unk']
        JOB_PREFIX: str = my_params['prefix_job']

    # if any error happens
    except Exception as exc:
        my_output_msg(format_exc())
        raise ValueError(
            "Incorrect or missing parameters file ; check path and/or contents of json reference.") from exc

    my_logs_global_config(LOG_FILE, verbose=args.verbose)

    core_call(args.exclude, args.multithreading, args.build, args.params,
              TAXAS_LEVELS, DATABASE, SAMPLES_PATH, JOB_PREFIX)
