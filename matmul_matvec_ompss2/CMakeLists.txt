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

if (Python_FOUND AND Python_Interpreter_FOUND AND Python_NumPy_FOUND)
  # install (copy) scripts
  configure_file(validate.py validate.py @ONLY)
else ()
  message("Python interpreter not found")
endif () # PYTHONINTERP_FOUND

# Function to create tests for executable.
function(create_test_pair executable)
  # Execution test
  add_test(NAME Test_${executable}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 256 16 1 1)
  # Create FIXTURES_SETUP for Test_${executable}
  set_tests_properties(Test_${executable} PROPERTIES
    ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml"
    FIXTURES_SETUP ST_${executable})

  # Add result execution test
  if (Python_FOUND AND Python_Interpreter_FOUND AND Python_NumPy_FOUND)
    # Add dependent test when Test_${executable} succeeded
    add_test(NAME Test_${executable}_results
      COMMAND ${Python_EXECUTABLE} validate.py
      ${executable}_A.mat ${executable}_B.mat ${executable}_C.mat)
    # Add dependencies to make the test run after Test_${executable}
    # and run only if that test succeded and created the output files.
    set_tests_properties(Test_${executable}_results PROPERTIES
      FIXTURES_REQUIRED ST_${executable}
      DEPENDS Test_${executable}
      REQUIRED_FILES ${executable}_A.mat ${executable}_B.mat ${executable}_C.mat)
  endif () # Python_FOUND
endfunction()

# Create tests
SET(TESTLIST matmul matvec) # The order is important
foreach (ISMATVEC 0;1)
  list(GET TESTLIST ${ISMATVEC} TEST)

  # weak tests.
  SET(FETCHLIST nofetch fetchall fetchfirst taskfor) # fetch
  foreach (FETCHTASK 0;1;2;3)
    list(GET FETCHLIST ${FETCHTASK} FETCHSUFIX)

    SET(NODELIST nonode node) # node clause
    foreach (WITH_NODE 0;1)
      list(GET NODELIST ${WITH_NODE} NODESUFIX)

      set(executable "${TEST}_weak_${FETCHSUFIX}_${NODESUFIX}")
      message("- Add benchmark: ${executable}")
      add_executable(${executable} matmul_ompss2.c)
      target_compile_definitions(${executable}
	PUBLIC -DISSTRONG=0
	PUBLIC -DFETCHTASK=${FETCHTASK}
	PUBLIC -DWITH_NODE=${WITH_NODE}
	PUBLIC -DISMATVEC=${ISMATVEC})
      target_link_libraries(${executable} argparser)

      create_test_pair(${executable})
    endforeach ()
  endforeach ()

  # Strong
  SET(STRONGLIST ignore flat nested)
  foreach (STRONGTYPE 1;2)
    list(GET STRONGLIST ${STRONGTYPE} STRONGSUFIX)

    # Executable
    set(executable "${TEST}_strong_${STRONGSUFIX}")
    message("- Add benchmark: ${executable}")
    add_executable(${executable} matmul_ompss2.c)
    target_compile_definitions(${executable}
      PUBLIC -DISSTRONG=${STRONGTYPE}
      PUBLIC -DISMATVEC=${ISMATVEC})
    target_link_libraries(${executable} argparser)

    create_test_pair(${executable})
  endforeach ()

  configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_${TEST}_dim.sh @ONLY)
  configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_${TEST}_dim.sh @ONLY)

endforeach ()

# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

# Process
configure_file("${PROJECT_SOURCE_DIR}/utils/process_mn.py" process_mn.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/process_dim.py" process_dim.py COPYONLY)

