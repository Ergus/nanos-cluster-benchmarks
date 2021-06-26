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

# 3.9 required to use FindOpenMP
# 3.12.4 required to use FindPython
cmake_minimum_required(VERSION 3.12.4)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(MERCURIUM REQUIRED COMPONENTS intel)

project(nanos6_benchmarks)

option(WITH_NAMESPACE "Build tests for namespace code" true)
option(WITH_ASAN "Build tests with address sanitizer" false)
option(WITH_EXTRAE "Build with Extrae traces support" false)
option(WITH_KEEP "Keep mcc intermediate files" false)
option(WITH_SERIAL "Build serial version of the code when available." false)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_compile_definitions(_GNU_SOURCE)
add_compile_options("--ompss-2" $<$<BOOL:WITH_KEEP>:-k>)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --ompss-2")

# Enable address sanitizer
if (WITH_ASAN)
  message("Using ASAN")
  add_compile_definitions(__WITH_ASAN=1)
  add_compile_options(-fno-omit-frame-pointer -fsanitize=address)

  set (CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
endif ()

# Enable Extrae
set(USE_INSTRUMENT "none")
if (WITH_EXTRAE)
  set(USE_INSTRUMENT "extrae")
  find_package(Extrae REQUIRED COMPONENTS nanosmpitrace)
  if (Extrae_FOUND)
    add_compile_definitions(__WITH_EXTRAE=1)
    include_directories(${EXTRAE_INCLUDES})
    link_libraries(${EXTRAE_LIBRARIES})
  endif ()
endif()

find_package(Numa REQUIRED)
include_directories(${NUMA_INCLUDES})
link_libraries(${NUMA_LIBRARIES})

find_package(MPI REQUIRED)
link_libraries(MPI::MPI_C)

find_package(Python MODULE COMPONENTS Interpreter NumPy)

find_package(INTEL REQUIRED)

# add the subdirectories here
add_subdirectory(utils/ArgParserC)

include_directories(utils)
enable_testing()

configure_file(${PROJECT_SOURCE_DIR}/utils/ArgParserBash/argparse.sh
  ${PROJECT_BINARY_DIR}/argparse.sh COPYONLY)

# Lists of include files
file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*_ompss2)
# Include all subdirectories with _ompss2 suffix.
foreach (child ${children})
  if ((IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
      AND (EXISTS ${CMAKE_CURRENT_LIST_DIR}/${child}/CMakeLists.txt))
    message ("-- Adding subdirectory: ${child}")
    add_subdirectory(${child})
  endif ()
endforeach ()

# For the toml file ==============
set(USE_NAMESPACE_TOML false)
if (WITH_NAMESPACE)
  message("-- Using NAMESPACE")
  set(USE_NAMESPACE_TOML true)
endif ()

set(VERSION_DEBUG_TOML true)
set(USE_TURBO_TOML false)
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(VERSION_DEBUG_TOML false)
  set(USE_TURBO_TOML true)
endif ()

configure_file(nanos6.toml nanos6.toml @ONLY)

