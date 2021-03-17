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

data = {}

# def load_file(input_file):
#     """Process the files and print the data in json format."""

#     # data[filename] = json.load(input_file)
#     print("Loading: ", filename)
#     data[filename] = pd.DataFrame(json.load(input_file))

#     return filename

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


def add_lines(ax1, ax2, dt_ts, column, label):
    assert dt_ts.empty == False
    dt_ts = dt_ts.sort_values(by=['worldsize'])

    x = dt_ts['worldsize']

    # Time graph
    y=dt_ts[column]
    err=dt_ts[column+' stdev'].divide(dt_ts['executions']**1/2)
    ax1.errorbar(x, y, err, fmt ='o-', label=label)

    fact = 1
    if dt_ts[x == fact].empty:
        fact = 2

    if dt_ts[x == fact].empty:
        print (" Can't create scalability graph for: ", label)
        return

    one = dt_ts[x == fact][column].values[0]

    # Scalability
    sy = one / y
    erry = sy * err / y
    ax2.errorbar(x, sy, erry, fmt ='o-', label=label)

def process_tasksize(rows, ts, column):
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
        prefix = key_prefix + " "

    fig.suptitle(prefix+str(rows) + " x " + str(ts))

    for key in data:
        dt_key = data[key]
        key_splits=key.split("_")

        if  key.endswith("mpi"):

            label = key_splits[-1]

            dt = dt_key[(dt_key['Rows'] == rows) & \
                        (dt_key['Tasksize'] == ts) ]

            add_lines_and_scale(ax1, ax2, dt, column, label)

        else:
            for ns in range(2):
                label = key_splits[1]+" "+[" ns", " nons"][ns]

                dt = dt_key[(dt_key['Rows'] == rows) & \
                            (dt_key['Tasksize'] == ts) & \
                            (dt_key['namespace_enabled'] == ns) ]

                add_lines_and_scale(ax1, ax2, dt, column, label)

    plt.legend(loc='upper right',
                   bbox_to_anchor=(1.14, 2.2),
                   fancybox=True, shadow=True, ncol=1)

    filename = column+"_"+str(rows)+"_"+str(ts)+".png"
    plt.savefig(filename)
    plt.close()
    print("Generated: ", filename)

def process_all_tasksize():
    """Create all the blocksize graphs"""

    first_dt = data[list(data.keys())[0]]

    rows_vals = first_dt['Rows'].drop_duplicates()
    ts_vals = first_dt['Tasksize'].drop_duplicates()

    for row in rows_vals:
        for ts in ts_vals:
            process_tasksize(row, ts, "Algorithm_time")


# for rows in row_vals:
#     dt_rows = dt[dt['Rows'] == rows]
#     tasksize_vals = dt_rows['Tasksize'].drop_duplicates()

#     fig = plt.figure()
#     gs = fig.add_gridspec(2, hspace=0)
#     (ax1, ax2) = gs.subplots(sharex=True, sharey=False)
#     fig.suptitle(key)

#     ax1.set_xlabel('Nodes')
#     ax1.set_ylabel(column)
#     ax1.grid(color='b', ls = '-.', lw = 0.25)

#     ax2.set_xlabel('Nodes')
#     ax2.set_ylabel("Scalability "+column)
#     ax2.grid(color='b', ls = '-.', lw = 0.25)

#     # Time
#     for ts in tasksize_vals:
#         dt_ts = dt_rows[dt_rows['Tasksize'] == ts]
#         if dt_ts.empty:
#             print("Ignore Tasksize: ", ts)
#             continue

#         add_lines(ax1, ax2, dt_ts, column, str(ts))


# def process_group(dim, ts, keys, prefix):
#     print("Processing group: ", dim, ts)
#     fig = plt.figure()
#     gs = fig.add_gridspec(2, hspace=0)
#     (ax1, ax2) = gs.subplots(sharex=True, sharey=False)
#     fig.suptitle(str(dim)+"x"+str(ts))

#     #First Graph (Time)
#     ax1.set_xlabel('Nodes')
#     ax1.set_ylabel("Algorithm time")
#     ax1.grid(color='b', ls = '-.', lw = 0.25)

#     #Second Graph (Scalability)
#     ax2.set_xlabel('Nodes')
#     ax2.set_ylabel("Scalability ")
#     ax2.grid(color='b', ls = '-.', lw = 0.25)

#     for key in keys:
#         if key.find("strong") != -1:
#             tag = "strong"
#         elif key.find("weak") != -1:
#             split = key.split("_")
#             tag = split[2]+" "+split[4]
#         elif key.find("mpi") != -1:
#             tag = "mpi"

#         dt = data[key]
#         dt_ts = dt.loc[(dt['Rows'] == dim) & (dt['Tasksize'] == ts)]

#         if dt_ts.empty:
#             print("Ignore key %s, Rows: %s, Tasksize: %s"
#                   % (key, dim, ts))
#             continue

#         add_lines(ax1, ax2, dt_ts, "Algorithm time", tag)

#     plt.legend(loc='center right',
#                bbox_to_anchor=(1.14, 2),
#                fancybox=True, shadow=True, ncol=1)

#     filename = prefix+"_"+str(dim)+"_"+str(ts)+".png"
#     plt.savefig(filename)
#     plt.close()
#     print("Generated: ", filename)


# def process_all(prefix, myregex):
#     dt = next(iter(data.values()))

#     row_vals = dt['Rows'].drop_duplicates()   # number of rows tests.
#     ts_vals = dt['Tasksize'].drop_duplicates()   # number of rows tests.

#     # Iterate over keys.
#     re1 = re.compile(myregex)    # Comments like # Anything
#     keys = [i for i in data if re1.search(i)]

#     for dim in row_vals:
#         for ts in ts_vals:
#             process_group(dim, ts, keys, prefix)
import pprint

if __name__ == "__main__":
    data = {}

    for fname in sys.argv[1:]:
        try:
            print("Loading:", fname)
            with open(fname, 'r') as f:
                fdata = json.load(f)

                for key in fdata:
                    df_in = pd.DataFrame(fdata[key])

                    if key in data:
                        data[key] = data[key].append(df_in)
                    else:
                        data[key] = df_in

        except IOError:
            print("File not accessible")

    process_all_tasksize()

    # process_all("Compare1", "mpi|(matvec_weak_fetchfirst)|(matvec_strong)")
    # process_all("Compare2", "(matvec_weak)|(matvec_strong)")

    # print(json.dumps(data, indent=4))

