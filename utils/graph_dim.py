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
import sys
import os
import json
import pandas as pd
import numpy as np
import re
import pprint
import pickle

def add_lines_and_scale(ax1, ax2, dt_ts, column, label):
    if dt_ts.empty:
        print("Ignoring: ", label)
        return

    assert dt_ts.empty == False
    dt_ts = dt_ts.sort_values(by=['worldsize'])
    x = dt_ts['worldsize']

    # Time graph
    y = dt_ts[column]
    erry=dt_ts[column+'_stdev'].divide(dt_ts['executions']**1/2)

    ax1.errorbar(x, y, erry, fmt ='o-', label=label)

    row_one = dt_ts[x == 1]
    if row_one.empty:
        print("No single node for : ", label)
        return

    one = row_one[column].values[0]
    errone = (row_one[column+'_stdev'] / (row_one['executions']**1/2)).values[0]

    # Scalability
    sy = one / dt_ts[column]
    errsy = sy * (erry/y + errone/one)

    ax2.errorbar(x, sy, errsy, fmt ='o-', label=label)


def process_tasksize(data, rows, ts, cpu_count, column):
    """Create graphs time vs tasksize"""

    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    (ax1, ax2) = gs.subplots(sharex=True, sharey=False)

    ax1.set_ylabel(column)
    ax1.grid(color='b', ls = '-.', lw = 0.25)

    ax2.set_xlabel('Nodes')
    ax2.set_ylabel("Scalability")
    ax2.grid(color='b', ls = '-.', lw = 0.25)

    prefix = ""
    keylist = list(data.keys())
    key_prefix = keylist[0].split("_")[0]

    if all(key.split("_")[0] == key_prefix for key in keylist):
        prefix = key_prefix

    fig.suptitle(prefix + " " + str(rows) + " x " + str(ts) + " (" + str(cpu_count) + " cores)")

    for key in data:
        dt_key = data[key]

        # "cholesky_memory_ompss2" -> "memory ompss2"
        label = " ".join(key.split("_")[1:])

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
    filename = column + "_" \
        + key.split("_")[0] + "_" \
        + str(rows) + "_" \
        + str(ts) + "_" \
        + str(cpu_count)


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

    plt.close()
    print("Generated: ", filename)


# Tasksize

def process_all(data):
    """Create all the blocksize graphs"""

    first_dt = data[list(data.keys())[0]]

    rows_vals = first_dt['Rows'].drop_duplicates()
    ts_vals = first_dt['Tasksize'].drop_duplicates()
    cpu_counts = first_dt['cpu_count'].drop_duplicates()

    for row in rows_vals:
        for ts in ts_vals:
            for cpu_count in cpu_counts:
                process_tasksize(data, row, ts, cpu_count, "Algorithm_time")


if __name__ == "__main__":
    data = {}

    if len(sys.argv) != 2:
        raise ValueError("Needs one filename as argument")

    fname = sys.argv[1]

    if fname.split(".")[-1] != "json":
        raise ValueError("Input file is not a .json")

    try:
        print("Loading:", fname)

        with open(fname, 'r') as f:
            fdata = json.load(f)

            for key in fdata:
                # key is the experiment name like: cholesky_fare_ompss2_taskfor

                df_in = pd.DataFrame(fdata[key])  # data to Pandas Dataframe

                if key in data:
                    data[key].append(df_in)
                else:
                    data[key] = df_in

    except IOError:
        print("File not accessible or json corrupt")

    process_all(data)

    # process_all("Compare1", "mpi|(matvec_weak_fetchfirst)|(matvec_strong)")
    # process_all("Compare2", "(matvec_weak)|(matvec_strong)")

    # print(json.dumps(data, indent=4))

