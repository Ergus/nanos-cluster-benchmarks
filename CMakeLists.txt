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

cmake_minimum_required (VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(MERCURIUM REQUIRED)

project(matvec)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(--ompss-2)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --ompss-2")

# add the subdirectories here
add_subdirectory(utils/ArgParserC)
include_directories(utils/ArgParserC)
include_directories(utils/cmacros)

configure_file(${PROJECT_SOURCE_DIR}/utils/ArgParserBash/argparse.sh
  ${PROJECT_BINARY_DIR}/argparse.sh COPYONLY)

file (GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)

foreach (child ${children})
  if ((IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
      AND (EXISTS ${CMAKE_CURRENT_LIST_DIR}/${child}/CMakeLists.txt))
    message ("-- Adding subdirectory: ${child}")
    add_subdirectory(${child})
  endif ()
endforeach ()
