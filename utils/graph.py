#!/usr/bin/env python3

# Copyright (C) 2021  Jimmy Aguilar Mena

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, os
import json
import pandas as pd
import numpy as np
import pprint
import pickle
from typing import *

import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.style'] = 'normal'
plt.rcParams['font.size'] = '8'


def save_all_files(filename: str, fig):
    """Save the graphs to two files."""

    # Save the plots with pickle to recover them:
    #
    # import matplotlib.pyplot as plt
    # import pickle
    # with open('filename.pkl', 'rb') as pkl:
    #     ax = pickle.load(pkl)
    # plt.show()
    with open(filename + ".pkl",'wb') as pkl:
        pickle.dump(fig, pkl)

    # Save the image as a png
    fig.savefig(filename + ".png",
                dpi=300,
                format='png',
                #bbox_extra_artists=[leg],
                bbox_inches='tight')

def add_lines(ax, dt_ts , column : str, label: str):
    "Add raw data graph"
    if dt_ts.empty:
        print("Ignoring: ", label)
        return

    dt_sorted = dt_ts.sort_values(by=['worldsize'])
    x = dt_sorted['worldsize']

    if dt_sorted[x == 1].empty:
        print("No single node for : ", label)
        return

    y = dt_sorted[column]
    erry = dt_sorted[column+'_stdev'].divide(dt_sorted['executions']**(1/2))

    ax.errorbar(x, y, yerr = erry, fmt = 'o-', linewidth=0.75, markersize=2, label=label)

def add_scalability(ax, dt_ts , column : str, label: str):
    """Add lines to the graphs."""
    if dt_ts.empty:
        print("Ignoring: ", label)
        return

    dt_sorted = dt_ts.sort_values(by=['worldsize'])
    x = dt_sorted['worldsize']

    row_one = dt_sorted[x == 1]
    if row_one.empty:
        print("No single node for : ", label)
        return

    one : float = row_one[column].values[0]
    errone : float = (row_one[column+'_stdev'] / (row_one['executions']**(1/2))).values[0]

    # Error
    y = dt_sorted[column]
    erry = dt_sorted[column+'_stdev'].divide(dt_sorted['executions']**(1/2))

    sy = one / dt_sorted[column]
    errsy = sy * (erry/y + errone/one)

    ax.errorbar(x, sy, errsy, fmt ='o-', linewidth=1, markersize=2, label=label)


def process_tasksize(data,
                     keyslist:list,
                     rows:int,
                     ts:int,
                     cpu_count:int,
                     column:str):
    """Create graphs time vs tasksize"""

    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    (ax1, ax2) = gs.subplots(sharex=True, sharey=False)

    ax1.set_ylabel(column)
    ax1.grid(color='b', ls = '-.', lw = 0.25)

    ax2.set_xlabel('Nodes')
    ax2.set_ylabel("Scalability")
    ax2.grid(color='b', ls = '-.', lw = 0.25)

    key_prefix : str = keyslist[0].split("_")[0]

    prefix: str = key_prefix if all(key.split("_")[0] == key_prefix for key in keyslist) else ""

    fig.suptitle(prefix + " " + str(rows) + " x " + str(ts) + " (" + str(cpu_count) + " cores)")

    for key in data:
        dt_key = data[key]
        # "cholesky_memory_ompss2" -> "memory ompss2"
        label : str = " ".join(key.split("_")[1:])

        dt = dt_key[(dt_key['Rows'] == rows) & \
                    (dt_key['Tasksize'] == ts) & \
                    (dt_key['cpu_count'] == cpu_count)]

        if key.endswith("mpi"):
            add_lines(ax1, dt, column, label)
            add_scalability(ax2, dt, column, label)
        else:
            for ns in range(2):
                dt_ns = dt[(dt['namespace_enabled'] == ns)]
                labelns = label + [" nons", " ns"][ns]

                add_lines(ax1, dt_ns, column, labelns)
                add_scalability(ax2, dt_ns, column, labelns)

    plt.legend(bbox_to_anchor=(1,1),
               loc='center left',
               fontsize='x-small',
               fancybox=True, shadow=True, ncol=1)

    # Save image file.
    filename = "Sacalability_" \
        + column + "_" \
        + key.split("_")[0] + "_" \
        + str(rows) + "_" \
        + str(ts) + "_" \
        + str(cpu_count)

    save_all_files(filename, fig)

    plt.close()
    print("Generated: ", filename)


def process_experiment(data, key:str,
                       rows:int,
                       ts_vals_list:list[int],
                       cpu_count_list:list[int],
                       column:str):
    """Create graphs comparing all the TS for same size and num_cpus"""

    fig = plt.figure()
    gs = fig.add_gridspec(nrows=3, ncols=(len(cpu_count_list)), hspace=0, wspace=0)
    ax = gs.subplots(sharex=True, sharey="row")

    # Title
    fig.suptitle(key + " " + str(rows))

    dt_key = data[key]

    nodes_list : list[int] = data[key]['worldsize'].drop_duplicates().sort_values().array

    for i in range(len(cpu_count_list)):

        cpu_count : int = cpu_count_list[i]
        ax1 = ax[0, i]
        ax2 = ax[1, i]
        ax3 = ax[2, i]

        if (i == 0):
            ax1.set_ylabel("Time")
            ax2.set_ylabel("Time (log)")
            ax3.set_ylabel("Scalability")

        ax1.title.set_text(str(cpu_count) + "cores")
        ax2.set_yscale('log')
        ax3.set_xlabel('Nodes')

        ax1.grid(color='b', ls = '-.', lw = 0.25)
        ax1.set_xticks(nodes_list)
        ax2.grid(color='b', ls = '-.', lw = 0.25)
        ax3.grid(color='b', ls = '-.', lw = 0.25)

        for ts in ts_vals_list:
            # "cholesky_memory_ompss2" -> "memory ompss2"
            label : str = str(ts)

            dt = dt_key[(dt_key['Rows'] == rows) & \
                        (dt_key['Tasksize'] == ts) & \
                        (dt_key['cpu_count'] == cpu_count)]

            if key.endswith("mpi"):
                add_lines(ax1, dt, column, label)
                add_lines(ax2, dt, column, label)
                add_scalability(ax3, dt, column, label)
            else:
                for ns in range(2):
                    dt_ns = dt[(dt['namespace_enabled'] == ns)]
                    labelns = label + [" nons", " ns"][ns]

                    add_lines(ax1, dt_ns, column, labelns)
                    add_lines(ax2, dt_ns, column, labelns)
                    add_scalability(ax3, dt_ns, column, labelns)

    plt.subplots_adjust()

    plt.legend(bbox_to_anchor=(1,1),
               loc='center left',
               fontsize='x-small',
               fancybox=True, shadow=True, ncol=1)

    # Save image file.
    # Save image file.
    filename = "Compare_" + column + "_" + key + "_" + str(rows)

    save_all_files(filename, fig)

    plt.close()
    print("Generated: ", filename)


# Tasksize
def process_all(data):
    """Create all the blocksize graphs"""

    keys_list : list[str] = list(data.keys());
    first_dt = data[keys_list[0]]

    # Get all the keys
    rows_list : list[int] = first_dt['Rows'].drop_duplicates().sort_values().array
    ts_list : list[int] = first_dt['Tasksize'].drop_duplicates().sort_values().array
    cpu_list : list[int] = first_dt['cpu_count'].drop_duplicates().sort_values().array

    # All benchmarks
    for row in rows_list:
        for ts in ts_list:
            for cpu_count in cpu_list:
                process_tasksize(data, keys_list, row, ts, cpu_count, "Algorithm_time")

    # All ts for same experiment.
    for key in keys_list:
        for row in rows_list:
            process_experiment(data, key, row, ts_list, cpu_list, "Algorithm_time")


if __name__ == "__main__":
    pd.set_option('display.max_rows', None)
    data : Dict[str, pd.DataFrame] = {}

    if len(sys.argv) < 2:
        raise ValueError("Needs one filename as argument")

    for fname in sys.argv[1:]:

        if fname.split(".")[-1] != "json":
            raise ValueError("Input file is not a .json")

        try:
            print("Loading:", fname)

            with open(fname, 'r') as f:
                fdata = json.load(f)

                for key in fdata:
                    # key is the experiment name like: cholesky_fare_ompss2_taskfor

                    df_in: pd.DataFrame = pd.DataFrame(fdata[key])  # data to Pandas Dataframe

                    if key in data:
                        data[key] = data[key].append(df_in)
                    else:
                        data[key] = df_in

        except IOError:
            print("File not accessible or json corrupt")

    process_all(data)

    # process_all("Compare1", "mpi|(matvec_weak_fetchfirst)|(matvec_strong)")
    # process_all("Compare2", "(matvec_weak)|(matvec_strong)")

    # print(json.dumps(data, indent=4))

