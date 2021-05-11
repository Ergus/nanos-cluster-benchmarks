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

###############################################################################
# Find Numa
#
# This sets the following variables:
# NUMA_FOUND        - True if NUMA found.
# NUMA_INCLUDES     - where to find numa.h, etc.
# NUMA_LIBRARIES    - List of libraries when using NUMA.
###############################################################################

include (FindPackageHandleStandardArgs)

set (NUMA_ROOT_DIR $ENV{NUMA_ROOT_DIR})
message (STATUS "NUMA_ROOT_DIR: ${NUMA_ROOT_DIR}")

# Find where are the header files.
find_path(NUMA_INCLUDE
  NAMES numa.h numaif.h
  HINTS ENV NUMA_ROOT_DIR
  PATHS ENV C_INCLUDE_PATH
  DOC "Numa include path"
  PATH_SUFFIXES include)

find_library(NUMA_LIBRARY
  NAMES ${Extrae_FIND_COMPONENTS}
  HINTS ENV NUMA_ROOT_DIR
  PATHS ENV LD_LIBRARY_PATH
  REQUIRED
  DOC "Numa library"
  PATH_SUFFIXES lib lib64 lib32)

find_package_handle_standard_args (NUMA
  REQUIRED_VARS NUMA_LIBRARY NUMA_INCLUDE)

mark_as_advanced(NUMA_INCLUDES NUMA_LIBRARIES)

if (NUMA_FOUND)
  set (NUMA_INCLUDES  ${NUMA_INCLUDE})
  set (NUMA_LIBRARIES ${NUMA_LIBRARY})
endif ()
