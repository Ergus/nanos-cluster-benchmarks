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

import sys
import os
import re
from statistics import mean, stdev

import json

re_ignore = re.compile('(# [^-=]+)|(^performance)')    # Comments like # Anything
re_pair = re.compile('(?P<key>\w+): (?P<value>.+)')   # KEY: value

re_number = re.compile('(?P<number>\d+(?P<float>\.\d+(e\+\d+)?)?)') # KEY: number

re_next = re.compile('# -+')            # Divisor like # ---------
re_report = re.compile('# =+')          # Divisor like # =========

results = {}

def process_group(a_dict):
    """Process group of executions with same parameters."""

    copydic = {}
    for key in a_dict:
        assert(isinstance(a_dict[key], list))
        vals_list = a_dict[key]
        count = len(vals_list)
        assert count > 0

        if "executions" in copydic:
            assert copydic["executions"] == count
        else:
            copydic["executions"] = count

        # keys int and string
        if all(isinstance(i, (str, int)) for i in vals_list) :
            if all(i == vals_list[0] for i in vals_list):
                copydic[key] = vals_list[0]
            else:
                print(vals_list)
                sys.exit("Failed `all` instances!")

        else:  # Floats
            if count == 1: # single element.
                assert isinstance(vals_list[0], float)
                m = vals_list[0]
                s = 0;
            else:          # Use mean and stdev when multiple elements
                assert all(isinstance(i, (float, int)) for i in vals_list)
                m = mean(vals_list)
                s = stdev(vals_list)

            copydic[key] = m
            copydic[key + "_stdev"] = s

    exe = copydic.pop("Executable")
    assert exe
    basename = os.path.basename(exe)
    results.setdefault(basename, []).append(copydic)

def process_file(input_file):
    """Process the files and print the data in json format."""

    line_dict = {}
    count = 0

    for line in input_file:
        # Ignore
        if re_ignore.match(line):
            continue

        # A pair value is always attached.
        # Even the filename for latter check
        pair = re_pair.match(line)
        if pair:
            key = pair.groupdict()['key']
            value = pair.groupdict()['value']

            number = re_number.match(value)
            if number:
                if number.groupdict()['float']: # it is a float so will be averaged later
                    value = float(number.groupdict()['number'])
                else:                           # it is a key so will be used as a key/info
                    value = int(number.groupdict()['number'])
            else:
                value = value.strip('\"')
            # Create or append
            line_dict.setdefault(key, []).append(value)
            continue

        # -------------- repetition
        if re_next.match(line):
            count = count + 1
            continue

        # ============== end group
        if re_report.match(line):
            if count > 0:
                process_group(line_dict)
                line_dict = {}
                count = 0
            continue

    if count > 0:
        process_group(line_dict)


if __name__ == "__main__":
    for dirname in sys.argv[1:]:
        if os.path.isdir(dirname):
            print("Going into:", dirname)
            results = {}
            # Process input
            for basename_in in os.listdir(dirname):
                if basename_in.endswith(".out"):
                    fname_in = os.path.join(dirname, basename_in)
                    print("Processing:", fname_in)
                    try:
                        with open(fname_in) as fin:
                            process_file(fin)
                    except IOError:
                        print("Couldn't read input:", fname_in)


            # Write output
            fname_out = dirname.rstrip(os.sep)+".json"
            print("Writing:", fname_out)
            try:
                with open(fname_out, "w") as fout:
                    json.dump(results, fout, indent=4)
            except IOError:
                print("Couldn't write output")

