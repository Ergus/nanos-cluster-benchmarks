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

SET(TESTTYPE weak strong)
foreach (TYPEID 0;1)
  list(GET TESTTYPE ${TYPEID} TYPE)
  set(executable "cholesky_fare_${TYPE}_ompss2")

  message("- Add benchmark: ${executable}")

  add_executable(${executable} cholesky_fare_ompss2.c)
  target_compile_definitions(${executable} PUBLIC -DTYPEID=${TYPEID})

  add_test(NAME Test_${executable}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 512 16)
  # Create FIXTURES_SETUP for Test_${executable}
  set_tests_properties(Test_${executable} PROPERTIES
    ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml")

endforeach ()

#======= cholesky_fare optimized versions ==========
SET(TESTTYPE memory taskfor)
foreach (TYPEID 0;1)
  list(GET TESTTYPE ${TYPEID} TYPE)
  set(executable "cholesky_fare_${TYPE}_ompss2")

  message("- Add benchmark: ${executable}")
  add_executable(${executable} "${executable}.c")

  add_test(NAME Test_${executable}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 512 16)
  set_tests_properties(Test_${executable} PROPERTIES
    ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml")

endforeach ()

# ==============================================
# submitter
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submitter_dim.sh" submitter_dim.sh @ONLY)

# Interactive
configure_file("${PROJECT_SOURCE_DIR}/utils/interactive_dim.sh" interactive_dim.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)
