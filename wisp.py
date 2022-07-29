#!/usr/bin/env python3.10
"This scripts aims to call for command line WISP over a set of files"
###########################################################################################
# This is your main interface. Or rather the call path to your main interface.
# Model parameters and such are stored inside a .json file, "wisp_params" by default
# If you want to change those, please refer to the README.md file!
###########################################################################################

from os import listdir, system
from sys import executable
from argparse import ArgumentParser
from traceback import format_exc
from shlex import split
from subprocess import Popen
from wisp_tools import my_output_msg, my_function_timer, my_logs_global_config, my_futures_collector
from wisp_lib import load_json


@my_function_timer("WISP run")
def core_call(method: str, exclusion: str, params: str, taxas_levels: list[str], db: str, unk_path: str, job_prefix: str, multithreading_state: int = 1) -> None:
    """
    Calls the building and prediction functions with global constants defined above
    If a job fails, skips to the next one
    """
    retcodes: list = []
    try:
        file_list: list[str] = listdir(unk_path)
        match method:
            case 'model':
                communicators = my_futures_collector(
                    func=Popen,
                    argslist=[[split(
                        f"{executable} wisp_model.py {db} {params} -l {[level]}")] for level in taxas_levels],
                    num_processes=min(multithreading_state,
                                      len(TAXAS_LEVELS)+1)
                )
            case 'predict':
                communicators = my_futures_collector(
                    func=Popen,
                    argslist=[[split(
                        f"{executable} wisp_predict.py {db} {params} {job_prefix}_{file[:-4]} -f {file} -e {exclusion}")] for file in file_list],
                    num_processes=multithreading_state
                )
            case 'build':
                communicators = my_futures_collector(
                    func=Popen,
                    argslist=[[split(
                        f"{executable} wisp_build.py {db} {params} -l {[level]}")] for level in taxas_levels],
                    num_processes=min(multithreading_state,
                                      len(TAXAS_LEVELS)+1)
                )
        retcodes = [p.communicate() for p in communicators]
    except Exception as excl:
        retcodes = ['Base error']
        raise BaseException("Job failed") from excl
    finally:
        for i, code in enumerate(retcodes):
            my_output_msg(f"Job {i} : {code}")
        system(f"rm -r data/{DATABASE}/temp")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "method", type=str, choices=['build', 'model', 'predict'], help="A callable func to execute")
    parser.add_argument(
        "params", help="path to a params .json file", type=str)
    parser.add_argument(
        "-v", "--verbose", help="Detailed log file", action='store_true')
    parser.add_argument(
        "-j", "--jobnumber", type=str, default='0', help="Gives a job identifier for WISP logs")
    parser.add_argument(
        "-e", "--exclude", type=str, help="Used for one-vs-all tests")
    parser.add_argument(
        "-l", "--leaveoneout", help="Leave the unknown sample out of training set", action='store_true')
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

    my_logs_global_config(
        filepath=LOG_FILE,
        identifier=args.jobnumber,
        verbose=args.verbose)

    parameters_path: str = f"parameters_files/{args.params}" if 'parameters_files' not in args.params else args.params

    core_call(
        method=args.method,
        exclusion=args.exclude,
        params=parameters_path if '.json' in parameters_path else f"{parameters_path}.json",
        taxas_levels=TAXAS_LEVELS,
        db=DATABASE,
        unk_path=SAMPLES_PATH,
        job_prefix=JOB_PREFIX
    )
