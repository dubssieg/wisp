from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler
import os
from pathlib import Path
import pickle
from typing import Any

import psutil


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
    mem_info = psutil.virtual_memory()
    load_avg = psutil.getloadavg()
    num_cpus = psutil.cpu_count()

    stats = {
        "load_avg_1min": load_avg[0],
        "load_avg_5min": load_avg[1],
        "load_avg_15min": load_avg[2],
        "num_cpus": num_cpus,
        "ram_usage_percent": mem_info.percent,
        "total_ram": mem_info.total,
        "used_ram": mem_info.used,
    }

    if as_str:
        return (
            f"Load Avg (1/5/15 min): {load_avg[0]:.2f}/{load_avg[1]:.2f}/{load_avg[2]:.2f}, "
            f"CPUs: {num_cpus}, "
            f"RAM Usage: {mem_info.percent:.2f}%, "
            f"Total RAM: {format_size(mem_info.total)}, "
            f"Used RAM: {format_size(mem_info.used)}"
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
        "%(asctime)s :: %(levelname)s :: %(name)s ::  %(process)d :: %(message)s"
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
    logger.addHandler(file_handler)

    for module in ignore_list:
        logging.getLogger(module).setLevel(logging.CRITICAL)


if __name__ == "__main__":
    print(system_stats())
    print(system_stats(as_str=True))
