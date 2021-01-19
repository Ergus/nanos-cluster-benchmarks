#!/usr/bin/python3

# Copyright (C) 2019  Jimmy Aguilar Mena

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

import numpy as np

A = np.loadtxt("matrix_A.mat")
x = np.loadtxt("vector_x.mat")
b = np.loadtxt("vector_b.mat")

b2 = np.dot(A, x)

norm = np.max(np.abs(b2 - b))

print ("Max difference %g" % norm)

if (norm > 1.0e-3):
    i = 0
    for xb, xb2 in zip(b, b2):
        if np.abs(xb - xb2) > 1.0e-3:
            print ("i = %d b[i] = %g b2[i] = %g" % (i, xb, xb2))
        i += 1
