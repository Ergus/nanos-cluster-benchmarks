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
# MERCURIUM_FOUND - True if MERCURIUM compiler was found in the path.
# MERCURIUM_CC  - Mercurium C compiler if found.
# MERCURIUM_CXX - Mercurium C++ compiler if found.
# MERCURIUM_FC  - Mercurium Fortran compiler if found.
###############################################################################

include(FindPackageHandleStandardArgs)

message(STATUS "CMAKE_MERCURIUM_PATH: ${CMAKE_MERCURIUM_PATH}")
message(STATUS "ENV_MERCURIUM_PATH: $ENV{MERCURIUM_PATH}")

unset(MERCURIUM_CC)
unset(MERCURIUM_CXX)
unset(MERCURIUM_FC)

# Search the different compilers
find_program(MERCURIUM_CC
  NAMES mcc
  HINTS ${CMAKE_MERCURIUM_PATH} ENV MERCURIUM_PATH
  DOC "MERCURIUM C compiler"
  PATH_SUFFIXES bin)

find_program(MERCURIUM_CXX
  NAMES mcxx
  HINTS ${CMAKE_MERCURIUM_PATH} ENV MERCURIUM_PATH
  DOC "MERCURIUM C++ compiler"
  PATH_SUFFIXES bin)

find_program(MERCURIUM_FC
  NAMES mfc
  HINTS ${CMAKE_MERCURIUM_PATH} ENV MERCURIUM_PATH
  DOC "MERCURIUM FORTRAN compiler"
  PATH_SUFFIXES bin)

find_package_handle_standard_args(MERCURIUM
  REQUIRED_VARS MERCURIUM_CC MERCURIUM_CXX MERCURIUM_FC)

mark_as_advanced(MERCURIUM_CC MERCURIUM_CXX MERCURIUM_FC)

if (MERCURIUM_FOUND)
  set(CMAKE_C_COMPILER ${MERCURIUM_CC})
  set(CMAKE_CXX_COMPILER ${MERCURIUM_CXX})
else ()
  message(WARNING "Mercurium was not found in default locations")
endif ()
