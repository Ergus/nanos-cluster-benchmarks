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

re_ignore = re.compile('# [^-=]+')    # Comments like # Anything
re_pair = re.compile('(?P<key>[\w,\s]+): (?P<value>\d+(?P<float>\.\d+(e\+\d+)?)?)') # KEY: number
re_next = re.compile('# -+')            # Divisor like # ---------
re_report = re.compile('# =+')          # Divisor like # =========
re_value = re.compile('((\w+)|performance$ time)$')

results = {}

def process_group(a_dict):
    """Process group of executions with same parameters."""
    for key in a_dict:
        if isinstance(a_dict[key], list):
            m = mean(a_dict[key])
            s = stdev(a_dict[key])
            a_dict[key] = (m, s)

    print(json.dumps(a_dict))


def process_file(input_file):
    '''Process the files and print the data in json format.'''
    line_dict = {}
    count = 0

    for line in input_file:
        # Ignore
        if re_ignore.match(line):
            continue

        # A pair value
        match = re_pair.match(line)
        if match:
            key = match.groupdict()['key']
            if match.groupdict()['float']: # it is a float so will be averaged later
                fvalue = float(match.groupdict()['value'])
                if key in line_dict:
                    line_dict[key].append(fvalue)
                else:
                    line_dict[key] = [fvalue]
            else:                  # it is a key so will be used as a key/info
                ivalue = int(match.groupdict()['value'])
                if key in line_dict:
                    assert line_dict[key] == ivalue
                    assert count > 0
                else:
                    line_dict[key] = ivalue

            continue

        # --------------
        if re_next.match(line):
            count = count + 1
            continue

        # ==============
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
        print ("Needs an input file")

