

from datetime import datetime

import psutil


def format_duration(seconds: float) -> str:
    """Convert seconds to a formatted duration string (hh:mm:ss)."""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def get_current_datetime_string() -> str:
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

def system_stats() -> dict:
    mem_info = psutil.virtual_memory()
    return {
        "cpu_usage" : psutil.cpu_percent(interval=1),        
        "ram_usage_percent" : mem_info.percent,
        "total_ram" : mem_info.total,
        "used_ram" : mem_info.used
    }
