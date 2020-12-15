#!/usr/bin/python3

import getopt, sys
import re
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

seq = {
    'program_name': './sequence_translator',
    'time_tokens': ['Reading time:', 'Computation time:', 'Writing time:', 'Total time:']
}

omp = {
    'program_name': './prog_OMP',
    'time_token': 'time ='
}

def plot_multi(data, cols=None, spacing=.1, **kwargs):

    from pandas import plotting

    # Get default color style from pandas - can be changed to any other color list
    if cols is None: cols = data.columns
    if len(cols) == 0: return
    colors = ['red', 'green', 'blue']

    # First axis
    ax = data.loc[:, cols[0]].plot(kind='scatter', label=cols[0], color=colors[0], **kwargs)
    ax.set_ylabel(ylabel=cols[0])
    lines, labels = ax.get_legend_handles_labels()

    for n in range(1, len(cols)):
        # Multiple y-axes
        ax_new = ax.twinx()
        ax_new.spines['right'].set_position(('axes', 1 + spacing * (n - 1)))
        data.loc[:, cols[n]].plot(kind='scatter', ax=ax_new, label=cols[n], color=colors[n % len(colors)], **kwargs)
        ax_new.set_ylabel(ylabel=cols[n])
        ax_new.set_ylim([0, data.loc[:, cols[n]].max()])

        # Proper legend position
        line, label = ax_new.get_legend_handles_labels()
        lines += line
        labels += label

    ax.legend(lines, labels, loc=0)
    ax.set_xlabel("Iteration")
    return ax

# [program] Dict
def collect_data(program, sizes):
    times = []
    thread_means = []
    for threads in range(1, 9):
        print("Threads:", threads)
        os.environ['OMP_NUM_THREADS'] = str(threads)
        for size in sizes:
            results = []
            for i in range(0, 3): 
                output = subprocess.run([program['program_name']] + str(size).split(), capture_output=True).stdout.decode('utf-8')
                result = re.search('(?<=' + program['time_token'] + ').\S*', output).group(0).strip()
                print((threads, size, result))
                results.append(float(result))
            results.remove(max(results))
            results.remove(min(results))
            mean = sum(results)/len(results)
            thread_means.append((threads, size, mean))
    return thread_means 

# [program] Dict
def collect_seq_data(program, sizes):
    times = []
    for i in range(0, 8): 
        output = subprocess.run([program['program_name']] + str(problem_size).split(), capture_output=True).stdout.decode('utf-8')
        times = [] # An array where each element is one time per execution
        for token in program['time_tokens']:
            result = re.search('(?<=' + token + ').\S*', output).group(0).strip()
            print(token, result)
            times.append(float(result))
        results.append(times)

    
    print("Sequential results:", 
    results = {}
    # zip(*) creates array of arrays where each array is a type of time, i.e. total time
    for time_type, time_set in zip(program['time_tokens'], zip(*results)): 
        print(time_type, ':', time_set)
        time_set.remove(max(time_set))
        time_set.remove(min(time_set))
        mean = sum(time_set)/len(time_set)
        results[time_type] = mean;

    return results


# index represents the index of sizes array
def calculate_speedup(baseline_mean, test_times, sizes):
    speedup_list = []
    # datum = (threads, problem_size, mean_time)
    print(baseline_mean, sizes)
    for datum in test_times:
        speedup = baseline_mean[2] / datum[2] * 100
        speedup_list.append((datum[0], datum[1], speedup))
    return speedup_list

# collect_seq_data
# returns a dict of key-value pairs
# each key is a value from program['time_types']
# value for each key is mean execution time
seq_means = collect_seq_data(seq)

omp_results = collect_data(omp)
omp_speedups = calculate_speedup(seq_mean, omp_results, sizes)


print(omp_speedups[0])
data = pd.DataFrame({
    "Threads": [i for i in range(1, 9)],
    "300": [speedup[2] for speedup in omp_speedups if speedup[1] == sizes[0]], 
    "3K": [speedup[2] for speedup in omp_speedups if speedup[1] == sizes[1]],
    "30K": [speedup[2] for speedup in omp_speedups if speedup[1] == sizes[2]]
})


filename = 'test'
data.to_csv(filename + '.csv', index = False)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.title('Speedup SEQ vs OMP')
ax.scatter(data['Threads'], data['300'], label = '300')
ax.scatter(data['Threads'], data['3K'], label = '3K')
ax.scatter(data['Threads'], data['30K'], label = '30K')
plt.xlabel('Threads')
plt.ylabel('Speedup (%)')
plt.legend()
plt.savefig('results.png')

plt.show()    

