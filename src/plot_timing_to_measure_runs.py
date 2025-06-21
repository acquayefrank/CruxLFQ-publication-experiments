import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CRUX_LOG_DIR = "../results/crux_runs_timing_memory_logs"
FLASH_LOG_DIR = "../results/flash2.0_runs_timing_memory_logs"
memory_usage_pattern = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")
time_pattern = re.compile(r"Elapsed \(wall clock\) time.*: (\d+:)?(\d+):([\d.]+)")

def parse_log_file(log_file):
    with open(log_file, 'r') as file:
        content = file.read()
    memory_usage = memory_usage_pattern.search(content)
    time_usage = time_pattern.search(content)
    # print(memory_usage, time_usage)
    if memory_usage and time_usage:
        memory = int(memory_usage.group(1)) / 1048576 #1024  # Convert to MB
        hours = int(time_usage.group(1)[:-1]) if time_usage.group(1) else 0
        minutes = int(time_usage.group(2))
        seconds = float(time_usage.group(3))
        time = (hours * 3600 + minutes * 60 + seconds)/60  # Convert to seconds
        return memory, time
    return None, None

def collect_data(log_dir):
    memory_data = {}
    time_data = {}
    for i in range(1, 4):
        for num_files in [1, 2, 4, 8, 16, 20]:
            log_file = os.path.join(log_dir, f"iteration_{i}_{num_files}.log")
            memory, time = parse_log_file(log_file)
            if memory is not None and time is not None:
                if num_files not in memory_data:
                    memory_data[num_files] = []
                    time_data[num_files] = []
                memory_data[num_files].append(memory)
                time_data[num_files].append(time)
    return memory_data, time_data

def plot_data(crux_memory_data, crux_time_data, flash_memory_data, flash_time_data):
    num_files = sorted(crux_memory_data.keys())
    crux_avg_memory = [sum(crux_memory_data[n]) / len(crux_memory_data[n]) for n in num_files]
    crux_avg_time = [sum(crux_time_data[n]) / len(crux_time_data[n]) for n in num_files]
    flash_avg_memory = [sum(flash_memory_data[n]) / len(flash_memory_data[n]) for n in num_files]
    flash_avg_time = [sum(flash_time_data[n]) / len(flash_time_data[n]) for n in num_files]
    x_ticks = range(0, 21, 2)

    plt.figure(figsize=(25, 17))

    plt.plot(num_files, flash_avg_memory, marker='o', label='FlashLFQ', markersize=15, linewidth=2)
    plt.plot(num_files, crux_avg_memory, marker='o', label='CruxLFQ', markersize=15, linewidth=2)
    
    plt.xlabel('Number of spectrum files', fontsize=55)
    plt.ylabel('Average memory usage (GB)', fontsize=55)
    # plt.title('Average Memory Usage vs Number of Files')
    plt.legend(fontsize=55, loc="best") # using a size in points
    plt.xticks(x_ticks, fontsize=55)  # Set the X-axis ticks to the specific number of files
    plt.yticks(fontsize=55)
    plt.tight_layout()
    plt.savefig('average_memory_usage.pdf')
    plt.close()

    plt.figure(figsize=(25, 17))
    plt.plot(num_files, flash_avg_time, marker='o', label='FlashLFQ', markersize=15, linewidth=2)
    plt.plot(num_files, crux_avg_time, marker='o', label='CruxLFQ', markersize=15, linewidth=2)
    
    # Linear regression for FlashLFQ time
    flash_time_slope, flash_time_intercept = np.polyfit(num_files, flash_avg_time, 1)
    flash_time_line = np.poly1d([flash_time_slope, flash_time_intercept])
    plt.plot(num_files, flash_time_line(num_files), label=f'FlashLFQ Fit: y={flash_time_slope:.2f}x+{flash_time_intercept:.2f}', linestyle='--')

    plt.xlabel('Number of spectrum files', fontsize=55)
    plt.ylabel('Average wall time (minutes)', fontsize=55)
    # plt.title('Average Time vs Number of Files')
    plt.legend(fontsize=55, loc="best") # using a size in points
    plt.xticks(x_ticks, fontsize=55)  # Set the X-axis ticks to the specific number of files
    plt.yticks(fontsize=55)
    plt.tight_layout()
    plt.savefig('average_time.pdf')
    plt.close()


def save_data_to_tsv(data, filename):
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(data, orient='index')
    df.index.name = 'Number of Files'
    df.columns = [f'Iteration {i+1}' for i in range(df.shape[1])]
    # Save the DataFrame to a TSV file
    df.to_csv(filename, sep='\t')


if __name__ == "__main__":
    crux_memory_data, crux_time_data = collect_data(CRUX_LOG_DIR)
    flash_memory_data, flash_time_data = collect_data(FLASH_LOG_DIR)

    # Save each dataset to a separate TSV file
    save_data_to_tsv(crux_memory_data, 'crux_memory_data.tsv')
    save_data_to_tsv(crux_time_data, 'crux_time_data.tsv')
    save_data_to_tsv(flash_memory_data, 'flash_memory_data.tsv')
    save_data_to_tsv(flash_time_data, 'flash_time_data.tsv')
    
    plot_data(crux_memory_data, crux_time_data, flash_memory_data, flash_time_data)