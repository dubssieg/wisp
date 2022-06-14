from time import monotonic
from datetime import timedelta
from logging import basicConfig, INFO, WARNING, info, warning
from typing import Callable


def my_output_msg(string: str, severity: Callable = info) -> None:
    "Date and time of action + info specified in str"
    severity(f"{string}")


def my_logs_clear(filepath: str):
    "Cleans out log file if run was success"
    with open(f"{filepath}", 'w'):
        pass


def my_logs_global_config(filepath: str, verbose: bool = False):
    "Clears log and defines logging info"
    my_logs_clear(f"{filepath}.log")
    basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', filename=f"{filepath}.log",
                encoding='utf-8', level=INFO if verbose else WARNING)


def my_function_timer(arg: str):
    """
    Decorator ; prints out execution time of decorated func
    * arg : descrptor name of job
    """
    def my_inner_dec(func):
        def wrapper(*args, **kwargs):
            my_output_msg("Starting job...")
            start_time = monotonic()
            res = func(*args, **kwargs)
            end_time = monotonic()
            my_output_msg(
                f"{arg} : Finished after {timedelta(seconds=end_time - start_time)} seconds")
            return res
        return wrapper
    return my_inner_dec
