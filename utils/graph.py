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

data = {}

def load_file(input_file):
    """Process the files and print the data in json format."""
    basename = os.path.basename(input_file.name)
    filename = os.path.splitext(basename)[0]

    # data[filename] = json.load(input_file)
    print("Loading: ", filename)
    data[filename] = pd.DataFrame(json.load(input_file))
    print(data[filename])

    # print(data[filename][data[filename]['Tasksize'] == 32])
    # print(data[filename]['Tasksize'].drop_duplicates())
    # row_vals = data[filename]['Rows'].drop_duplicates().to_list()

    # print(data[filename]['Rows'].drop_duplicates().to_list())

    return filename

def process_tasksize(key, column):
    """Create graphs time vs tasksize"""
    dt = data[key]
    row_vals = dt['Rows'].drop_duplicates()   # number of rows tests.

    for rows in row_vals:
        dt_rows = dt[dt['Rows'] == rows]
        tasksize_vals = dt_rows['Tasksize'].drop_duplicates()

        fig = plt.figure()
        gs = fig.add_gridspec(2, hspace=0)
        (ax1, ax2) = gs.subplots(sharex=True, sharey=False)
        fig.suptitle(key)

        ax1.set_xlabel('Nodes')
        ax1.set_ylabel(column)
        ax1.grid(color='b', ls = '-.', lw = 0.25)

        ax2.set_xlabel('Nodes')
        ax2.set_ylabel("Scalability "+column)
        ax2.grid(color='b', ls = '-.', lw = 0.25)

        # Time
        for ts in tasksize_vals:

            dt_ts = dt_rows[dt_rows['Tasksize'] == ts]

            row_one = dt_ts[dt_ts['worldsize'] == 1]
            one = row_one[column].values[0] / row_one['Iterations'].values[0]

            # Time graph
            x=dt_ts['worldsize']
            y=dt_ts[column].divide(dt_ts['Iterations'])
            err=dt_ts[column+' stdev'].divide(dt_ts['Iterations']*dt_ts['executions']**1/2)
            ax1.errorbar(x, y, err, fmt ='o-', label=str(ts))

            # Scalability
            sy = one / y
            erry = sy * err / y
            ax2.errorbar(x, sy, erry, fmt ='o-', label=str(ts))

        plt.legend(loc='center right',
                   bbox_to_anchor=(1.14, 1),
                   fancybox=True, shadow=True, ncol=1)

        filename = column.replace(" ", "_")+"_"+key+"_"+str(rows)+".png"
        plt.savefig(filename)
        plt.close()
        print("Generated: ", filename)

def process_bars_tasksize(key):
    """Create graphs time vs tasksize"""
    dt = data[key]
    ts_vals = dt['Tasksize'].drop_duplicates()   # number of rows tests.

    for ts in ts_vals:
        dt_ts = dt[dt['Tasksize'] == ts]   # rows with tasksize == ts

        row_vals = dt_ts['Rows'].drop_duplicates()

        # Figure
        fig, ax = plt.subplots()

        fig.suptitle("Tasksize " + str(ts))

        ax.set_xlabel('Nodes')
        ax.set_ylabel('Time')
        ax.grid(axis = 'y', color = 'b', ls = '-.', lw = 0.25)

        nbars = len(row_vals)
        width = 1.0 / (nbars + 1)
        offset = - width * nbars / 4 # /2 center /2 half

        for rows in row_vals:
            dt_rows = dt_ts[dt_ts['Rows'] == rows]

            x = np.arange(dt_rows['worldsize'].size)

            ax.bar(x + offset, dt_rows["Total time"], color = 'b',  width = width)
            ax.bar(x + offset, dt_rows["Algorithm time"], color = 'r', width = width)
            offset = offset + width

        x = dt_ts[dt_ts['Rows'] == rows]['worldsize']
        ax.set_xticks(np.arange(x.size))
        ax.set_xticklabels([str(i) for i in x])

        filename = "Bars_"+key+"_"+str(ts)+".png"
        plt.savefig(filename)
        plt.close()
        print("Generated: ", filename)


if __name__ == "__main__":
    for fname in sys.argv[1:]:
        try:
            with open(fname, 'r') as f:
                key = load_file(f)
                process_tasksize(key, "Algorithm time")
                process_tasksize(key, "Total time")
                process_bars_tasksize(key)

        except IOError:
            print("File not accessible")

    # print(json.dumps(data, indent=4))

