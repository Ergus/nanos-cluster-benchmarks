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
SET(NODELIST nonode node)
foreach (NODEID 0;1)
  list(GET NODELIST ${NODEID} NODESUFIX)

  if ((NOT WITH_NONODE) AND (NODEID EQUAL 0))
    message("-- Skipping no-node")
    continue()
  endif()

  SET(TESTLIST matmul matvec) # The order is important
  foreach (ISMATVEC 0;1)
    list(GET TESTLIST ${ISMATVEC} TEST)

    set(BLASLIST noblas blas)
    foreach(BLAS 0 1)
      list(GET BLASLIST ${BLAS} BLASSUFIX)

      # Strong
      SET(TASKTYPELIST strong_flat_task strong_nested_task weak_taskfor)
      foreach (TASKTYPEID 0;1;2)
	list(GET TASKTYPELIST ${TASKTYPEID} TASKTYPESUFIX)

	# Executable
	set(executable "${TEST}_${TASKTYPESUFIX}_${NODESUFIX}_${BLASSUFIX}_ompss2")
	message("- Add benchmark: ${executable}")
	add_executable(${executable} matmul_ompss2.c)
	target_compile_definitions(${executable}
	  PUBLIC WITHNODE=${NODEID}
	  PUBLIC TASKTYPE=${TASKTYPEID}
	  PUBLIC ISMATVEC=${ISMATVEC}
	  PUBLIC BLAS=${BLAS})
	target_link_libraries(${executable} argparser)

	create_test_pair(${executable})
      endforeach ()

      # weak tests.
      SET(FETCHLIST weak_nofetch_task weak_fetchall_task) # fetch
      foreach (FETCHID 0;1)
	list(GET FETCHLIST ${FETCHID} FETCHSUFIX)

	set(executable "${TEST}_${FETCHSUFIX}_${NODESUFIX}_${BLASSUFIX}_ompss2")
	message("- Add benchmark: ${executable}")
	add_executable(${executable} matmul_ompss2.c)
	target_compile_definitions(${executable}
	  PUBLIC TASKTYPE=3
	  PUBLIC WITHNODE=${NODEID}
	  PUBLIC ISMATVEC=${ISMATVEC}
	  PUBLIC FETCHTASK=${FETCHID} 
	  PUBLIC BLAS=${BLAS})

	create_test_pair(${executable})
      endforeach ()
    endforeach ()
  endforeach ()
endforeach ()

# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_dim.sh @ONLY)

# Interactive
configure_file("${PROJECT_SOURCE_DIR}/utils/interactive_dim.sh" interactive_dim.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

