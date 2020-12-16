#!/usr/bin/python3

import getopt, sys
import re
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fasta_file = 'allChromosomes_pretty.fa'

seq = {
    'program_name': './sequence_translator',
    'time_tokens': ['Reading time:', 'Computation time:', 'Writing time:', 'Total time:']
}

# OPT for Optimized, uses both MPI and OMP
opt = {
    'program_name': 'sequence_translatorOPT',
    'time_tokens': ['Reading time:', 'Computation time:', 'Writing time:', 'Total time:']
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
def collect_data(program, cores, multi=False):
    thread_mean_dicts = []

    exec_args = [program['program_name'], fasta_file, '../out.fa']
    # Launch across multiple machines if multi
    if multi:
        exec_args = ['mpirun', '-np', str(cores), '--hostfile', 'h2', '--mca', 'btl_tcp_if_include', 'eno1'] + exec_args
    else:
        exec_args = ['mpirun', '-np', str(cores)] + exec_args

    
    for threads in range(1, 9):
        print("Threads:", threads)
        os.environ['OMP_NUM_THREADS'] = str(threads)
        
        results = []
        for i in range(8): 
            output = subprocess.run(exec_args, capture_output=True).stdout.decode('utf-8')

            times = [] # Each elem is one time per exec
            for token in program['time_tokens']:
                result = re.search('(?<=' + token + ').\S*', output).group(0).strip()
                print(threads, token, result)
                times.append(float(result))
            results.append(times)
        
        print("Results for", threads, "threads:")
        results = {}
        # zip(*) creates array of arrays where each array is a type of time, i.e. total time
        for time_type, time_set in zip(program['time_tokens'], zip(*results)):
            print(time_type, ':', time_set)
            time_set.remove(max(time_set))
            time_set.remove(min(time_set))
            mean = sum(time_set)/len(time_set)
            results[time_type] = mean; 

        thread_mean_dicts.append((threads, results))

    return thread_mean_dicts 

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


# Returns [(threads, {'time_token': speedup, ...}), ...]
def calculate_speedup(baseline_means, test_times):
    speedups = []
    print(baseline_means, test_times)

    # datum = (threads, {token_string: mean_time, ...})
    for datum in test_times:
        speedup_dict = {}
        for key, val in datum[1]:
            speedup = baseline_means[key] / val * 100
            speedup_dict[key] = speedup

        speedups.append((datum[0], speedup_dict))
    return speedups

# collect_seq_data
# returns {'time_token': mean_time, ...}
# i.e., a dict of key-value pairs
# each key is a value from program['time_types']
# value for each key is mean execution time
seq_means = collect_seq_data(seq)

# Execute on one machine with 6 cores
opt6_results = collect_data(opt, 6, False)
# Execute using host file and 12 cores
opt12_results = collect_data(opt, 12, True)

# List of tuples where each tuple is (threads, {'time_token': speedup, ...})
opt6_speedups = calculate_speedup(seq_means, opt6_results)
opt12_speedups = calculate_speedup(seq_means, opt12_results)

print(opt6_speedups[0])
opt6_data = pd.DataFrame({
    "Threads": [i for i in range(1, 9)],
    "Reading Time": [tup['Reading time:'] for tup in opt6_speedups], 
    "Computation Time": [tup['Computation time:'] for tup in opt6_speedups], 
    "Writing Time": [tup['Writing time:'] for tup in opt6_speedups], 
    "Total Time": [tup['Total time:'] for tup in opt6_speedups] 
})


opt12_data = pd.DataFrame({
    "Threads": [i for i in range(1, 9)],
    "Reading Time": [tup['Reading time:'] for tup in opt12_speedups], 
    "Computation Time": [tup['Computation time:'] for tup in opt12_speedups], 
    "Writing Time": [tup['Writing time:'] for tup in opt12_speedups], 
    "Total Time": [tup['Total time:'] for tup in opt12_speedups] 
})


filename = 'Translator'
data.to_csv(filename + '.csv', index = False)

# ['Reading time:', 'Computation time:', 'Writing time:', 'Total time:']


fig = plt.figure()
ax = fig.add_subplot(111)
plt.title('Speedup SEQ vs OPT6')
ax.scatter(opt6_data['Threads'], opt6_data['Reading Time'], label = 'Reading Time')
ax.scatter(opt6_data['Threads'], opt6_data['Computation Time'], label = 'Computation Time')
ax.scatter(opt6_data['Threads'], opt6_data['Writing Time'], label = 'Writing Time')

ax.scatter(opt6_data['Threads'], opt6_data['Total Time'], label = 'Total Time')
plt.xlabel('Threads')
plt.ylabel('Speedup (%)')
plt.legend()
plt.savefig('opt6-results.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.title('Speedup SEQ vs OPT12')
ax.scatter(opt12_data['Threads'], opt12_data['Reading Time'], label = 'Reading Time')
ax.scatter(opt12_data['Threads'], opt12_data['Computation Time'], label = 'Computation Time')
ax.scatter(opt12_data['Threads'], opt12_data['Writing Time'], label = 'Writing Time')

ax.scatter(opt12_data['Threads'], opt12_data['Total Time'], label = 'Total Time')
plt.xlabel('Threads')
plt.ylabel('Speedup (%)')
plt.legend()
plt.savefig('opt12-results.png')

plt.show()    

