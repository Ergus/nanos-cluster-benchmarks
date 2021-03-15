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

results = []

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
            assert all(i == vals_list[0] for i in vals_list)
            copydic[key] = vals_list[0]

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

    results.append(copydic)

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
            key = match.groupdict()['key']
            value = match.groupdict()['value']

            number = re_pair.match(value)
            if number:
                if number.groupdict()['float']: # it is a float so will be averaged later
                    value = float(number.groupdict()['number'])
                else:                           # it is a key so will be used as a key/info
                    value = int(match.groupdict()['number'])

            if key in line_dict:                # append or create
                line_dict[key].append(fvalue)
            else:
                line_dict[key] = [fvalue]

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
    for fname in sys.argv[1:]:
        try:
            with open(fname) as f:
                # basename = os.path.splitext(sys.argv[1])[0]
                process_file(f)

        except IOError:
            print("File not accessible")
    else:
        print(json.dumps(results, indent=1))

