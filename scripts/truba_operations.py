import logging, psutil, csv


def print_memory_usage(verbose=False):
    
    if verbose:
        # Get the current process
        process = psutil.Process()

        # Get the memory info for the current process
        mem_info = process.memory_info()

        # Get the total memory and used memory in bytes
        total_memory = psutil.virtual_memory().total
        used_memory = mem_info.rss

        # Calculate the memory usage percentage
        memory_percentage = (used_memory / total_memory) * 100

        # Log the memory usage information
        logging.info(f"Used Memory: {used_memory / 1024 / 1024:.2f} MB, memory percentage: {memory_percentage:.2f}%")
        
        


def save_results_to_disk(results, result_file, verbose=False):
    """
    Save the results to a file.

    Args:
        results (list): A list of tuples, where each tuple contains the result of a job and the binary value.
        result_file (str): The path to the file where the results should be saved.
        verbose (bool): If True, log the message when the results are saved.
    """
    with open(result_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(results)

    if verbose:
        logging.info(f"Results saved to {result_file}")
