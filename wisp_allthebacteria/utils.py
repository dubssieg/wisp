from datetime import datetime

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


def system_stats() -> dict:
    mem_info = psutil.virtual_memory()
    return {
        "cpu_usage": psutil.cpu_percent(interval=1),
        "ram_usage_percent": mem_info.percent,
        "total_ram": mem_info.total,
        "used_ram": mem_info.used,
    }


if __name__ == "__main__":
    print(format_duration(100))
    print(format_duration(1000))
    print(format_duration(10000))
    print(format_duration(100000))
    print(format_duration(1000000))
    print(format_size(100))
    print(format_size(10000))
    print(format_size(1000000000))
    print(format_size(100000000000000))
    print(format_size(100000000000000000))
