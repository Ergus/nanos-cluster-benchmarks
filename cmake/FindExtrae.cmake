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
# Find Extrae
#
# This sets the following variables:
# EXTRAE_FOUND        - True if Extrae was found.
# EXTRAE_INCLUDES     - Directories containing the Extrae include files.
# EXTRAE_LIBRARIES    - Libraries needed to use Extrae.
###############################################################################

include (FindPackageHandleStandardArgs)

set (EXTRAE_HOME $ENV{EXTRAE_HOME})
message (STATUS "EXTRAE_HOME: ${EXTRAE_HOME}")

if (EXTRAE_HOME)

  # Find where are the header files.
  find_path (EXTRAE_INCLUDE
    NAMES extrae_user_events.h
    HINTS ENV EXTRAE_HOME
    PATHS ENV C_INCLUDE_PATH
    DOC "Extrae include path"
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)

  find_library (EXTRAE_LIBRARY
    NAMES ${Extrae_FIND_COMPONENTS}
    HINTS ENV EXTRAE_HOME
    PATHS ENV LD_LIBRARY_PATH
    REQUIRED
    DOC "Extrae library"
    PATH_SUFFIXES lib lib64 lib32
    NO_DEFAULT_PATH)

  find_package_handle_standard_args (Extrae
    REQUIRED_VARS EXTRAE_LIBRARY EXTRAE_INCLUDE)

endif ()

mark_as_advanced(EXTRAE_INCLUDES EXTRAE_LIBRARIES)

if (EXTRAE_FOUND)
  set(EXTRAE_INCLUDES  ${EXTRAE_INCLUDE})
  set(EXTRAE_LIBRARIES ${EXTRAE_LIBRARY})
endif ()
