#!/usr/bin/python3

import getopt, sys
import re
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fasta_file = '../allChromosomes_pretty.fa'

seq = {
    'program_name': './sequence-translator',
    'time_tokens': ['Final reading time:', 'Final computation time:', 'Final writing time:', 'Final total time:']
}

# OPT for Optimized, uses both MPI and OMP
opt = {
    'program_name': 'sequence-translatorOPT',
    'time_tokens': ['Final reading time:', 'Final computation time:', 'Final writing time:', 'Final total time:']
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

    
    #for threads in range(1, 9):
    #    print("Threads:", threads)
    #    os.environ['OMP_NUM_THREADS'] = str(threads)
        
    results = []
    for i in range(3): 
        print(exec_args)
        output = subprocess.run(exec_args, capture_output=True).stdout.decode('utf-8')
        print("OUTPUT", output)

        times = [] # Each elem is one time per exec
        for token in program['time_tokens']:
            print('Searching for', '(<=' + token + ').\S*')
            result = re.search('(?<=' + token + ').\S*', output).group(0).strip()
            print(token, result)
            times.append(float(result))
        results.append(times)
     
    print(results)   
#    print("Results for", threads, "threads:")
    # zip(*) creates array of arrays where each array is a type of time, i.e. total time
    print("zip results:", list(zip(*results)))
    means_dict = {}
    for time_type, time_list  in zip(program['time_tokens'], list(zip(*results))):
        time_list = list(time_list)
        print(time_type, ':', time_list)
        time_list.remove(max(time_list))
        time_list.remove(min(time_list))
        mean = sum(time_list)/len(time_list)
        means_dict[time_type] = mean; 

    print("Final means dict", means_dict)
    return means_dict 

# [program] Dict
def collect_seq_data(program):
    results = []
    for i in range(3): 
        output = subprocess.run([program['program_name'], '../allChromosomes_pretty.fa', '../out.fa'], capture_output=True).stdout.decode('utf-8')
        times = [] # An array where each element is one time per execution
        for token in program['time_tokens']:
            print('OUTPUT:', output)
            print('Searching for', '(<=' + token + ').\S*')
            result = re.search('(?<=' + token + ').\S*', output).group(0).strip()
            print(token, result, 'actual result:', result[:-1])
            times.append(float(result[:-1]))
        results.append(times)
    print("Seq times:", results)

    times_dict = {}
    # zip(*) creates array of arrays where each array is a type of time, i.e. total time
    for time_type, time_list in zip(program['time_tokens'], zip(*results)): 
        time_list = list(time_list)
        print(time_type, ':', time_list)
        time_list.remove(max(time_list))
        time_list.remove(min(time_list))
        mean = sum(time_list)/len(time_list)
        times_dict[time_type] = mean;

    print("Final seq results", times_dict)
    return times_dict


# Returns [(threads, {'time_token': speedup, ...}), ...]
def calculate_speedup(baseline_means, test_times):
    speedup_dict = {}
    print(baseline_means, test_times)

    # datum = (threads, {token_string: mean_time, ...})
    print("test_times:", test_times)
    for key, val in test_times.items():
        speedup = baseline_means[key] / val * 100
        speedup_dict[key] = speedup

    print("Final speedups:", speedup_dict)
    return speedup_dict

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
"""
opt6_data = pd.DataFrame({
    "Reading Time": opt6_speedups['Final reading time:'],
    "Computation Time": opt6_speedups['Final computation time:'], 
    "Writing Time": opt6_speedups['Final writing time:'], 
    "Total Time": opt6_speedups['Final total time:']
})


opt12_data = pd.DataFrame({
    "Reading Time": opt12_speedups['Final reading time:'],
    "Computation Time": opt12_speedups['Final computation time:'], 
    "Writing Time": opt12_speedups['Final writing time:'], 
    "Total Time": opt12_speedups['Final total time:'] 
})
"""

filename = 'Translator'

# ['Reading time:', 'Computation time:', 'Writing time:', 'Total time:']


fig = plt.figure()
plt.title('OPT6 vs OPT12 Speedups')
ax = fig.add_subplot(111)
X = np.arange(4)
opt6 = [value for _,value in opt6_speedups.items()] 
opt12 = [value for _,value in opt12_speedups.items()] 
p1 = plt.bar(X + 0.00, opt6, color = 'b', width = 0.25)
p2 = plt.bar(X + 0.25, opt12, color = 'g', width = 0.25)

 
#ax.scatter(opt6_, opt6_data['Reading Time'], label = 'Reading Time')
#ax.scatter(opt6_data['Threads'], opt6_data['Computation Time'], label = 'Computation Time')
#ax.scatter(opt6_data['Threads'], opt6_data['Writing Time'], label = 'Writing Time')

#ax.scatter(opt6_data['Threads'], opt6_data['Total Time'], label = 'Total Time')
#plt.xlabel('Threads')
plt.ylabel('Speedup (%)')
plt.xlabel('Time Metric')
plt.xticks(X, ('Reading Time', 'Computation Time', 'Writing Time', 'Total Time'))
plt.legend((p1[0], p2[0]), ('OPT6', 'OPT12'))
plt.savefig('opt-results.png')

#fig = plt.figure()
#ax = fig.add_subplot(111)
#plt.title('Speedup SEQ vs OPT12')
#ax.scatter(opt12_data['Threads'], opt12_data['Reading Time'], label = 'Reading Time')
#ax.scatter(opt12_data['Threads'], opt12_data['Computation Time'], label = 'Computation Time')
#ax.scatter(opt12_data['Threads'], opt12_data['Writing Time'], label = 'Writing Time')

#ax.scatter(opt12_data['Threads'], opt12_data['Total Time'], label = 'Total Time')
#plt.xlabel('Threads')
#plt.ylabel('Speedup (%)')
#plt.legend()
#plt.savefig('opt12-results.png')

plt.show()    

