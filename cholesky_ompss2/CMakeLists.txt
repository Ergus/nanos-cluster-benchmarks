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

set(BLA_VENDOR Intel10_64ilp)
find_package(BLAS REQUIRED)

find_path(
  MKL_INCLUDE_DIR mkl.h
  HINTS ENV MKLROOT
  PATH_SUFFIXES include)

# Add the executables and link
include_directories(${MKL_INCLUDE_DIR})

add_executable(cholesky_weak cholesky_weak.c cholesky_utils.c)

target_link_options(cholesky_weak PUBLIC ${BLAS_LINKER_FLAGS})

target_link_libraries(cholesky_weak argparser ${BLAS_LIBRARIES})

# # install (copy) scripts
# configure_file(submit_mn.sh submit_mn.sh @ONLY)
# configure_file(submiter_mn.sh submiter_mn.sh @ONLY)
# configure_file(trace.sh trace.sh @ONLY)
# configure_file(extrae_template.xml extrae.xml @ONLY)

# configure_file(process.py process.py COPYONLY)
# configure_file(create_paraver_trace.py create_paraver_trace.py COPYONLY)

# configure_file(validate.py validate.py @ONLY)
