#!/usr/bin/env python3

# Copyright (C) 2022  Jimmy Aguilar Mena

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

import os, sys, re

text :str = ""

try:
    with open(sys.argv[1], 'r') as f:
        text = f.read()

except IOError:
    print("File not accessible")

# Pattern for function
func :str = r'^void +(\w+)_ompss2\s*\([^)]*\)\s*' + r"({(?s:.+?)^})"
re_func = re.compile(func, flags=re.MULTILINE)

# Extract printf (where is the benchmark name.)
extract :str = r'^ *printf *\((.+?)\);'
re_extract = re.compile(extract, flags=re.MULTILINE)

matches = re.findall(re_func, text)

for i in matches:
    untab :str = re.sub(r'\t', ' ', i[1])
    noempty :str = re.sub(r'\n *\n', '\n', untab)
    nobreak :str = re.sub(r' *\\\n *', ' ', noempty)
    nocomments1 :str = re.sub(r' *//.*\n', '\n', nobreak) # //
    nocomments2 :str = re.sub(r' */\*(?s:.+?)\*/', "", nocomments1) # /* */

    print(nocomments)
    printmatch = re.search(re_extract, nocomments2)

    print(printmatch.group(0), printmatch.group(1))
