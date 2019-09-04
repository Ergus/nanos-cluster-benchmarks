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

###############################################################################
# Find libErep
#
# This sets the following variables:
# LIBEREP_FOUND        - True if liberep was found.
# LIBEREP_INCLUDES     - Directories containing the liberep include files.
# LIBEREP_LIBRARIES    - Libraries needed to use liberep.
###############################################################################

INCLUDE(FindPackageHandleStandardArgs)

MESSAGE(STATUS "CMAKE_LIBEREP_ROOT: ${CMAKE_LIBEREP_ROOT}")
MESSAGE(STATUS "LIBEREP_ROOT: $ENV{LIBEREP_ROOT}")

# Set extended paths
SET(INC_PATHS
  /usr/include
  /usr/local/include
  /usr/include )

SET(LIB_PATHS
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/lib/x86_64-linux-gnu)

# Find where are the header files.
FIND_PATH(LIBEREP_INCLUDE
  NAMES CommandLine.hpp
  HINTS ${CMAKE_LIBEREP_ROOT} ENV LIBEREP_ROOT
  PATHS ${INC_PATHS} ENV CPLUS_INCLUDE_PATH
  DOC "LIBEREP include path"
  PATH_SUFFIXES erep
  )

FIND_LIBRARY(LIBEREP_LIBRARY
  NAMES erep
  HINTS ${CMAKE_LIBEREP_ROOT} ENV LIBEREP_ROOT
  PATHS ${LIB_PATHS} ENV LD_LIBRARY_PATH
  DOC "LIBEREP library"
  PATH_SUFFIXES lib lib64 lib32)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBEREP
  REQUIRED_VARS LIBEREP_INCLUDE LIBEREP_LIBRARY)

MARK_AS_ADVANCED(LIBEREP_INCLUDE LIBEREP_LIBRARY)

