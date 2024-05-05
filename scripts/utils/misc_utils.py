import psutil
import time
import os
import json
import pandas as pd

from scripts.config import OUTPUT_DIR

def time_it(func):


    def wrapper(*args, **kwargs):
        process = psutil.Process()

        if hasattr(os, 'getloadavg'):
            load_before = os.getloadavg()  # System load average before function execution
        else:
            load_before = (None, None, None)

        start_time = time.time()  # Start timing
        memory_before = process.memory_info().rss / (1024 ** 2)

        result = func(*args, **kwargs)  # Call the function

        end_time = time.time()  # End timing
        memory_after = process.memory_info().rss / (1024 ** 2)
        # Calculate duration in minutes
        duration = (end_time - start_time) / 60
        cpu_core_used = process.cpu_num()

        if hasattr(os, 'getloadavg'):
            load_after = os.getloadavg()  # System load average after function execution
        else:
            load_after = (None, None, None)

        # Prepare a concise representation of arguments for JSON output
        args_repr = [f"{len(a)}" if isinstance(
            a, pd.DataFrame) else repr(a) for a in args]

        # Create a dictionary to store function name, duration, and arguments
        log_entry = {
            "function_name": func.__name__,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(start_time)),
            "duration_minutes": duration,
            "arguments": args_repr,
            "memory_before": memory_before,
            "memory_after": memory_after,
            "cpu_core_used": cpu_core_used,
            "load_before": load_before,
            "load_after": load_after
        }

        # Write to JSON file
        json_file = os.path.join(OUTPUT_DIR, 'function_timings.json')
        with open(json_file, 'a') as f:
            json.dump(log_entry, f)
            f.write('\n')  # Write a newline to separate entries

        return result
    return wrapper
