# Copyright (C) 2020  Jimmy Aguilar Mena

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


###############################################################################
# Find Mercurium
#
# This sets the following variables:
# INTEL_FOUND - True if Intel compiler was found in the path.
# ICC         - Mercurium C compiler if found.
# ICPC        - Mercurium C++ compiler if found.
# IFORT       - Mercurium Fortran compiler if found.
###############################################################################

include(FindPackageHandleStandardArgs)

message(STATUS "CMAKE_INTEL_PATH: ${CMAKE_INTEL_PATH}")
message(STATUS "ENV_INTEL_PATH: $ENV{INTEL_PATH}")

# Search the different compilers
find_program(INTEL_CC
  NAMES icc
  HINTS ${CMAKE_INTEL_PATH} ENV INTEL_PATH
  DOC "Intel C compiler"
  PATH_SUFFIXES bin)

find_program(INTEL_ICPC
  NAMES icpc
  HINTS ${CMAKE_INTEL_PATH} ENV INTEL_PATH
  DOC "Intel C++ compiler"
  PATH_SUFFIXES bin)

find_program(INTEL_FC
  NAMES ifort
  HINTS ${CMAKE_INTEL_PATH} ENV INTEL_PATH
  DOC "Intel FORTRAN compiler"
  PATH_SUFFIXES bin)

find_package_handle_standard_args(INTEL
  REQUIRED_VARS INTEL_CC INTEL_ICPC INTEL_FC)

mark_as_advanced(INTEL_CC INTEL_ICPC INTEL_FC)

if (INTEL_FOUND)
  set(ICC ${INTEL_CC})
  set(ICPC ${INTEL_ICPC})
  set(IFORT ${INTEL_FC})
else ()
  message(WARNING "Mercurium was not found in default locations")
endif ()
