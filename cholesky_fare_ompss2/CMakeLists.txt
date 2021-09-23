# Compilation


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

set(TEST "cholesky")

set(BLA_VENDOR Intel10_64lp_seq)
find_package(BLAS REQUIRED)

find_path(
  MKL_INCLUDE_DIR mkl.h
  HINTS ENV MKLROOT
  PATH_SUFFIXES include)

# Add the executables and link
include_directories(${MKL_INCLUDE_DIR})

SET(TESTTYPE weak strong)
foreach (TYPEID 0;1)
  list(GET TESTTYPE ${TYPEID} TYPE)
  set(executable "cholesky_fare_${TYPE}_ompss2")

  message("- Add benchmark: ${executable}")

  add_executable(${executable} cholesky_fare_ompss2.c)
  target_compile_definitions(${executable} PUBLIC -DTYPEID=${TYPEID})
  target_link_options(${executable} PUBLIC ${BLAS_LINKER_FLAGS})
  target_link_libraries(${executable} ${BLAS_LIBRARIES} argparser)

  add_test(NAME Test_${executable}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 512 16)
  # Create FIXTURES_SETUP for Test_${executable}
  set_tests_properties(Test_${executable} PROPERTIES
    ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml")

endforeach ()

# Taskfor
set(executable "cholesky_fare_taskfor_ompss2")
message("- Add benchmark: ${executable}")
add_executable(${executable} cholesky_fare_ompss2_taskfor.c)
target_link_options(${executable} PUBLIC ${BLAS_LINKER_FLAGS})
target_link_libraries(${executable} ${BLAS_LIBRARIES} argparser)

add_test(NAME Test_${executable}
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 512 16)
# Create FIXTURES_SETUP for Test_${executable}
set_tests_properties(Test_${executable} PROPERTIES
ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml")

# submiter
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_${TEST}_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_${TEST}_dim.sh @ONLY)

# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

# Process
configure_file("${PROJECT_SOURCE_DIR}/utils/process_mn.py" process_mn.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/process_dim.py" process_dim.py COPYONLY)

