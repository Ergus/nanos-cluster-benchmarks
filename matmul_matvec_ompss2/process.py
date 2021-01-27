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

re_ignore = re.compile('# [^-=]+')    # Comments like # Anything
re_pair=re.compile('([\w,\s]+): (\d+(\.\d+(e\+\d+)?)?)') # KEY: number
re_next=re.compile('# -+')            # Divisor like # ---------
re_report=re.compile('# =+')                # Divisor like # =========


def process_file(input_file):
    '''Process the files and print the data in json format.'''
    line_dict = {}
    count = 0

    for line in input_file:
        # Ignore
        if re_ignore.match(line):
            continue

        # A pair
        match = re_pair.match(line)
        if match:
            if match.groups()[2]:
                line_dict[match.groups()[0]] = float(match.groups()[1])
            else:
                line_dict[match.groups()[0]] = int(match.groups()[1])
            continue

        # --------------
        if re_next.match(line):
            count = count + 1
            print (line_dict)
            line_dict.clear()
            continue

        # ==============
        if re_report.match(line):
            if count > 0:
                count = 0

if __name__ == "__main__":
    if (len(sys.argv) >= 1):
        try:
            with open(sys.argv[1]) as f:
                # basename = os.path.splitext(sys.argv[1])[0]
                process_file(f)

        except IOError:
            print("File not accessible")
    else:
        print ("Needs an input file")

