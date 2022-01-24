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

import matplotlib.pyplot as plt
import sys, os
import json
import pandas as pd
import numpy as np
import pprint
import pickle
from typing import *


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


def add_lines_and_scale(ax1, ax2, dt_ts , column : str, label: str):
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

    # Time graph
    y = dt_sorted[column]
    erry = dt_sorted[column+'_stdev'].divide(dt_sorted['executions']**(1/2))

    ax1.errorbar(x, y, erry, fmt ='o-', label=label)

    one : float = row_one[column].values[0]
    errone : float = (row_one[column+'_stdev'] / (row_one['executions']**(1/2))).values[0]

    # Scalability
    sy = one / dt_sorted[column]
    errsy = sy * (erry/y + errone/one)

    ax2.errorbar(x, sy, errsy, fmt ='o-', label=label)


def process_tasksize(data, keyslist:list, rows:int, ts:int, cpu_count:int, column:str):
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
            add_lines_and_scale(ax1, ax2, dt, column, label)
        else:
            for ns in range(2):
                dt_ns = dt[(dt['namespace_enabled'] == ns)]
                labelns = label + [" nons", " ns"][ns]

                add_lines_and_scale(ax1, ax2, dt_ns, column, labelns)

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
                       cpu_count:int,
                       column:str):
    """Create graphs comparing all the TS for same size and num_cpus"""

    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    (ax1, ax2) = gs.subplots(sharex=True, sharey=False)

    ax1.set_ylabel(column)
    ax1.grid(color='b', ls = '-.', lw = 0.25)

    ax2.set_xlabel('Nodes')
    ax2.set_ylabel("Scalability")
    ax2.grid(color='b', ls = '-.', lw = 0.25)

    # Title
    fig.suptitle(key + " " + str(rows) + " (" + str(cpu_count) + " cores)")

    dt_key = data[key]

    for ts in ts_vals_list:
        # "cholesky_memory_ompss2" -> "memory ompss2"
        label : str = str(ts)

        dt = dt_key[(dt_key['Rows'] == rows) & \
                    (dt_key['Tasksize'] == ts) & \
                    (dt_key['cpu_count'] == cpu_count)]

        if key.endswith("mpi"):
            dt_mpi = dt
            add_lines_and_scale(ax1, ax2, dt_mpi, column, label)
        else:
            for ns in range(2):
                dt_ns = dt[(dt['namespace_enabled'] == ns)]
                labelns = label + [" nons", " ns"][ns]

                add_lines_and_scale(ax1, ax2, dt_ns, column, labelns)

    plt.legend(bbox_to_anchor=(1,1),
               loc='center left',
               fontsize='x-small',
               fancybox=True, shadow=True, ncol=1)

    # Save image file.
    # Save image file.
    filename = "Compare_" \
        + column + "_" \
        + key + "_" \
        + str(rows) + "_" \
        + str(cpu_count)

    save_all_files(filename, fig)

    plt.close()
    print("Generated: ", filename)


# Tasksize
def process_all(data):
    """Create all the blocksize graphs"""

    keys : list = list(data.keys());
    first_dt = data[keys[0]]

    # Get all the keys
    rows_vals : list = first_dt['Rows'].drop_duplicates()
    ts_vals : list = first_dt['Tasksize'].drop_duplicates()
    cpu_counts : list = first_dt['cpu_count'].drop_duplicates()

    # All benchmarks
    for row in rows_vals:
        for ts in ts_vals:
            for cpu_count in cpu_counts:
                process_tasksize(data, keys, row, ts, cpu_count, "Algorithm_time")

    # All ts for same experiment.
    for key in keys:
        for row in rows_vals:
            for cpu_count in cpu_counts:
                process_experiment(data, key, row, ts_vals, cpu_count, "Algorithm_time")


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

