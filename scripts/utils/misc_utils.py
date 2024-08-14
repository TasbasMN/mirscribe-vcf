import psutil
import time
import os
import json
import pandas as pd

from scripts.config import OUTPUT_DIR

def time_it(func=None, enabled=True):
    def decorator(f):
        if not enabled:
            return f

        def wrapper(*args, **kwargs):
            process = psutil.Process()

            if hasattr(os, 'getloadavg'):
                load_before = os.getloadavg()
            else:
                load_before = (None, None, None)

            start_time = time.time()
            memory_before = process.memory_info().rss / (1024 ** 2)

            result = f(*args, **kwargs)

            end_time = time.time()
            memory_after = process.memory_info().rss / (1024 ** 2)
            duration = (end_time - start_time) / 60
            cpu_core_used = process.cpu_num()

            if hasattr(os, 'getloadavg'):
                load_after = os.getloadavg()
            else:
                load_after = (None, None, None)

            args_repr = [f"{len(a)}" if isinstance(a, pd.DataFrame) else repr(a) for a in args]

            log_entry = {
                "function_name": f.__name__,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(start_time)),
                "duration_minutes": duration,
                "arguments": args_repr,
                "memory_before": memory_before,
                "memory_after": memory_after,
                "cpu_core_used": cpu_core_used,
                "load_before": load_before,
                "load_after": load_after
            }

            json_file = os.path.join(OUTPUT_DIR, 'function_timings.json')
            with open(json_file, 'a') as file:
                json.dump(log_entry, file)
                file.write('\n')

            return result
        return wrapper

    if func:
        return decorator(func)
    return decorator
