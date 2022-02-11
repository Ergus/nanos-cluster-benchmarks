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

# Add the executables and link

set(TEST "jacobi")

# CMAKE_C_FLAGS_RELEASE is a space separated string.
# But add_custom_command needs a list. In cmake lists are strings with ;
string(REGEX REPLACE "[ \t]+" ";" LIBBASE_FLAGS ${CMAKE_C_FLAGS_RELEASE})

# libbase.a
add_custom_command(
  OUTPUT libbase.o
  COMMAND ${ICC} -c ${CMAKE_CURRENT_SOURCE_DIR}/jacobi_base.c ${LIBBASE_FLAGS} -o libbase.o
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/jacobi_base.c
  COMMENT "Building base object")

add_library(base STATIC libbase.o)

SET(FETCHLIST nofetch fetchall)
foreach (FETCHID 0;1)
  list(GET FETCHLIST ${FETCHID} FETCHSUFIX)

  set(EXECUTABLE "${TEST}_task_${FETCHSUFIX}_ompss2")
  message("- Add benchmark: ${EXECUTABLE}")

  add_executable(${EXECUTABLE} jacobi_task_ompss2.c)
  target_compile_definitions(${EXECUTABLE} PUBLIC FETCHTASK=${FETCHID})
  target_link_libraries(${EXECUTABLE} base)

  add_test(NAME Test_${EXECUTABLE}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${EXECUTABLE} 128 8 3)
  set_tests_properties(Test_${EXECUTABLE} PROPERTIES
    ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml")
endforeach ()

# jacobi taskfor
set(EXECUTABLE "${TEST}_taskfor_ompss2")
message("- Add benchmark: ${EXECUTABLE}")

add_executable(${EXECUTABLE} jacobi_taskfor_ompss2.c)
target_link_libraries(${EXECUTABLE} argparser base)

add_test(NAME Test_${EXECUTABLE}
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${EXECUTABLE} 128 8 3)

# jacobi serial
if (WITH_SERIAL)
  set(EXECUTABLE "${TEST}_serial")
  message("- Add benchmark: ${EXECUTABLE}")

  add_executable(${EXECUTABLE} jacobi_serial.c)
  target_link_libraries(${EXECUTABLE} argparser base)

  add_test(NAME Test_jacobi_serial
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${EXECUTABLE} 128 8 3)
endif ()

# install (copy) scripts
# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_dim.sh @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

#Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)
