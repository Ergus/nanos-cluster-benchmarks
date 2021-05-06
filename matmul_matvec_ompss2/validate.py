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
import numpy as np

if __name__ == "__main__":

    A = np.loadtxt(sys.argv[1])
    print("A = ", A.shape)
    B = np.loadtxt(sys.argv[2])
    print("B = ", B.shape)
    C1 = np.loadtxt(sys.argv[3])
    print("C = ", C1.shape)

    C2 = np.dot(A, B)

    norm = np.max(np.abs(C2 - C1))

    if (norm > 1.0e-3):
        print("Error")
        i = 0
        for xb, xb2 in zip(C1, C2):
            if np.abs(xb - xb2) > 1.0e-3:
                print ("i = %d b[i] = %g b2[i] = %g" % (i, xb, xb2))
                i += 1
    else:
        print("Ok")

    print("Max difference: %g" % norm)

