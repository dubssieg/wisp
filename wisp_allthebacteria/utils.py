import concurrent.futures
import hashlib
import logging
import os
import pickle
import threading
import zlib
from datetime import datetime
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Any, Literal

import lz4.frame
import msgpack
import numpy as np
import psutil
import pympler
import pympler.asizeof
from tqdm.auto import tqdm

LOG = logging.getLogger(__name__)


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


def system_stats(pid: int = None, as_str: bool = False) -> dict:
    if pid is None:
        pid = os.getpid()
    try:
        current_process = psutil.Process(pid)
    except psutil.NoSuchProcess:
        LOG.warning(f"Process with PID {pid} no longer exists.")
        return {} if not as_str else "Process no longer exists."

    # Current process and children RAM
    try:
        mem_info = current_process.memory_full_info()
    except psutil.NoSuchProcess:
        LOG.warning(f"Process with PID {pid} no longer exists.")
        return {} if not as_str else "Process no longer exists."

    script_memory_usage = mem_info.uss
    for child in current_process.children(recursive=True):
        try:
            script_memory_usage += child.memory_full_info().uss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

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


def compress(
    data: Any, as_list: bool = False, format: Literal["zlib", "msgpack"] = "msgpack"
) -> bytes | list[bytes]:
    if as_list:
        return [compress(i, as_list=False, format=format) for i in data]
    if format == "zlib":
        return zlib.compress(pickle.dumps(data))
    elif format == "msgpack":
        return lz4.frame.compress(msgpack.packb(data, use_bin_type=True))
    else:
        raise ValueError(format)


def decompress(
    data: bytes | list[bytes],
    as_list: bool = False,
    format: Literal["zlib", "msgpack"] = "msgpack",
) -> Any:
    if as_list:
        return [decompress(i, as_list=False, format=format) for i in data]
    if format == "zlib":
        return pickle.loads(zlib.decompress(data))
    elif format == "msgpack":
        return msgpack.unpackb(lz4.frame.decompress(data), raw=False)
    else:
        raise ValueError(format)


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
    if not logger.hasHandlers():
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            # "%(asctime)s :: %(levelname)s :: %(name)s.%(funcName)s[%(lineno)s] :: %(process)d :: %(message)s"
            "%(asctime)s :: %(levelname)s :: %(name)s.%(funcName)s[%(lineno)s] :: %(message)s"
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


class FunctionLogger:
    def __init__(self, interval: float, func, *args, level=logging.DEBUG, **kwargs):
        """Log the result of a function call at regular intervals."""
        self._interval = interval
        self._func = func
        self._args = args
        self._kwargs = kwargs
        self._level = level
        self._stop_event = threading.Event()
        self._thread = threading.Thread(target=self._log_function_result, daemon=True)

    def _log_function_result(self):
        while not self._stop_event.is_set():
            result = self._func(*self._args, **self._kwargs)
            LOG.log(self._level, result)
            self._stop_event.wait(self._interval)

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._stop_event.set()
        self._thread.join()


class SystemStatsLogger(FunctionLogger):
    def __init__(
        self, interval: float = 60.0, pid: int | None = None, level=logging.DEBUG
    ):
        """Log CPU/RAM usage for current script using FunctionLogger."""
        super().__init__(interval, system_stats, pid, as_str=True, level=level)


def cleanup_zombie_processes(pid: int | None = None):
    """Kill all zombie processes remaining"""

    if pid is None:
        pid = psutil.Process().pid
    for child in psutil.Process(pid).children(recursive=True):
        if child.status() == psutil.STATUS_ZOMBIE:
            LOG.warning(f"Killing zombie process {child.pid}")
            try:
                child.terminate()
                child.wait(5)  # Attendre 5 secondes pour qu'il se ferme
            except psutil.NoSuchProcess:
                pass
            except psutil.TimeoutExpired:
                LOG.error(
                    f"Failed to terminate zombie process {child.pid}, forcing kill"
                )
                child.kill()


def hash(data: Any) -> str:
    if not isinstance(data, (tuple, list, np.array)):
        data = [data]
    hasher = hashlib.sha256()
    for item in sorted(map(str, data)):
        hasher.update(item.encode())
    return hasher.hexdigest()


def get_weights(counts: list[int], balance_factor: float) -> list[float]:
    """Balancing sample generator.
    - balancing_factor=0.0 -> original distribution"
    - balancing_factor=1.0 -> balanced distribution"
    """
    if not counts:
        return []
    total_samples = sum(counts)
    if balance_factor == 0:
        return [count / total_samples for count in counts]
    else:
        uniform_weight = 1 / len(counts)
        return [
            (1 - balance_factor) * (count / total_samples)
            + balance_factor * uniform_weight
            for count in counts
        ]


def sample_count_estimation(counts: list[int], balance_factor: float) -> int:
    "On average, how many samples can we expect?"
    weights = get_weights(counts, balance_factor)
    weighted_samples = [
        count / weight if weight > 0 else float("inf")
        for count, weight in zip(counts, weights)
    ]
    return int(min(weighted_samples))


def sizeof(obj: Any, detail: bool = False) -> int | str:
    if detail:
        return pympler.asizeof.asized(obj, detail=1).format()
    return pympler.asizeof.asizeof(obj)
