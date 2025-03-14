from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
import os
from pathlib import Path
import pickle
from typing import Any
import concurrent.futures
import psutil
from tqdm.auto import tqdm


def format_duration(seconds: float) -> str:
    days = int(seconds // (24 * 3600))
    seconds %= 24 * 3600
    hours = int(seconds // 3600)
    seconds %= 3600
    minutes = int(seconds // 60)
    seconds = int(seconds % 60)

    hms = f"{hours:02}:{minutes:02}:{seconds:02}"
    if days > 0:
        return f"{days} jours + {hms}"
    else:
        return hms


def format_size(size_bytes: int) -> str:
    units = ["bytes", "KiB", "MiB", "GiB", "TiB"]
    unit_size = 1024

    for unit in units[:-1]:
        if size_bytes < unit_size:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= unit_size

    return f"{size_bytes:.2f} {units[-1]}"


def get_current_datetime_string() -> str:
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


def system_stats(as_str: bool = False) -> dict:
    current_process = psutil.Process(os.getpid())

    # Current process and children RAM
    mem_info = current_process.memory_full_info()
    script_memory_usage = sum(
        (
            child.memory_full_info().uss
            for child in current_process.children(recursive=True)
        ),
        start=mem_info.uss,
    )

    # System total RAM
    total_system_memory = psutil.virtual_memory().total

    # % Script RAM
    script_memory_percent = (script_memory_usage / total_system_memory) * 100

    # CPUs
    load_avg = psutil.getloadavg()
    num_cpus = psutil.cpu_count()

    stats = {
        "load_avg_1min": load_avg[0],
        "load_avg_5min": load_avg[1],
        "load_avg_15min": load_avg[2],
        "num_cpus": num_cpus,
        "script_memory_usage": script_memory_usage,
        "script_memory_percent": script_memory_percent,
        "total_system_memory": total_system_memory,
    }

    if as_str:
        return (
            f"RAM (Script): {format_size(script_memory_usage)} ({script_memory_percent:.2f}%), "
            f"CPU (1/5/15 min): {load_avg[0]:.2f}/{load_avg[1]:.2f}/{load_avg[2]:.2f}"
            # f"CPUs: {num_cpus}, "
            # f"Total RAM: {format_size(total_system_memory)}"
        )

    return stats


def serialize(obj: Any, path: str | Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as file:
        pickle.dump(obj, file)


def deserialize(path: str | Path) -> Any:
    path = Path(path)
    with path.open("rb") as file:
        return pickle.load(file)


def cpu_count(needed: str | int = 8) -> int:
    if isinstance(needed, int):
        return needed
    max_cpu = os.cpu_count()
    if max_cpu is None:
        raise RuntimeError("Cannot get cpu count")
    if needed == "max":
        return max_cpu
    elif needed == "max_minus_1" and max_cpu > 1:
        return max_cpu - 1
    elif needed == "max_minus_2" and max_cpu > 2:
        return max_cpu - 2
    elif needed == "half":
        return max_cpu // 2
    else:
        raise RuntimeError(f"Cannot get cpu count with: {needed}")


class BrokenProcessPoolFilter(logging.Filter):
    def filter(self, record):
        if record.exc_info:
            exc_type = record.exc_info[0]
            if issubclass(exc_type, concurrent.futures.process.BrokenProcessPool):
                record.exc_text = None
        return True


def config_logger(
    log_path: str | Path,
    terminal_level: str,
    file_level: str,
    file_size: int,
    file_count: int,
    ignore_list: list[str],
) -> None:
    """Terminal + files configuration."""
    # common config
    logger = logging.getLogger("")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s :: %(levelname)s :: %(name)s.%(funcName)s[%(lineno)s] :: %(process)d :: %(message)s"
    )

    # terminal config
    terminal_handler = logging.StreamHandler()
    terminal_handler.setFormatter(formatter)
    terminal_handler.setLevel(getattr(logging, terminal_level.upper()))
    logger.addHandler(terminal_handler)

    # files config
    path = Path(log_path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = RotatingFileHandler(
        path,
        mode="a",
        maxBytes=file_size,
        backupCount=file_count,
    )
    file_handler.setFormatter(formatter)
    file_handler.setLevel(getattr(logging, file_level.upper()))
    file_handler.addFilter(BrokenProcessPoolFilter())
    logger.addHandler(file_handler)

    for module in ignore_list:
        logging.getLogger(module).setLevel(logging.CRITICAL)


def slurm_tqdm(iterable, *args, **kwargs):
    # usage slurm_tqdm(..., disable=True)
    if os.getenv("SLURM_JOB_ID") and kwargs.pop("disable", False):
        return iterable
    return tqdm(iterable, *args, **kwargs)


def space_format(number: int):
    return f"{number:_}".replace("_", " ")


if __name__ == "__main__":
    print(system_stats())
    print(system_stats(as_str=True))
