from time import monotonic
from datetime import timedelta, datetime
from enum import Enum
from logging import basicConfig, INFO, info
from inspect import signature


def my_output_msg(string: str) -> None:
    "Prints the date and time of action + info specified in str"
    info(f"{string}")
    print(f"[{str(datetime.now())}] {string}")


def my_logs_clear(filepath: str):
    "Cleans out log file if run was success"
    with open(f"{filepath}", 'w'):
        pass


def my_logs_global_config(filepath: str):
    basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', filename=f"{filepath}.log",
                encoding='utf-8', level=INFO)


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


class Unk(Enum):
    pass


class my_Enum(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name


def my_entries_checker(func):
    """
    Permet de check les états des paramètres donnés en entrée à une fonction
    """
    def wrapper(*args, **kwargs):
        """
        Corps de test
        """

        # formation des listes de comparaison
        arguments = args + tuple(kwargs.keys())
        args_types = [type(arg) for arg in arguments]
        signature_func = signature(func).parameters
        annotations = [
            signature_func[elt].annotation for elt in signature_func]
        my_output_msg(f"Function {func.__name__} :")
        my_output_msg(f"Input params : {annotations}")
        # execution de la fonction décorée
        retour = func(*args, **kwargs)

        return retour
    return wrapper
