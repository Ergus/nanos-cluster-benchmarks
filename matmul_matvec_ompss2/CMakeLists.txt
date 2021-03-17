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

# The order is important
SET(TESTLIST matmul matvec)
foreach (ISMATVEC 0;1)
  list(GET TESTLIST ${ISMATVEC} TEST)

  # Weak
  SET(FETCHLIST nofetch fetchall fetchfirst)
  foreach (FETCHTASK 0;1;2)
    list(GET FETCHLIST ${FETCHTASK} FETCHSUFIX)

    SET(NODELIST nonode node)
    foreach (WITH_NODE 0;1)
      list(GET NODELIST ${WITH_NODE} NODESUFIX)

      set(executable "${TEST}_weak_${FETCHSUFIX}_${NODESUFIX}")
      message("Add benchmark: ${executable}")
      add_executable(${executable} matmul_ompss2.c)
      target_compile_definitions(${executable}
	PUBLIC -DFETCHTASK=${FETCHTASK}
	PUBLIC -DWITH_NODE=${WITH_NODE}
	PUBLIC -DISMATVEC=${ISMATVEC})
      target_link_libraries(${executable} argparser)

    endforeach ()
  endforeach ()

  # Strong
  SET(STRONGLIST ignore flat interleave nested)
  foreach (STRONGTYPE 1;2;3)
    list(GET STRONGLIST ${STRONGTYPE} STRONGSUFIX)

    set(executable "${TEST}_strong_${STRONGSUFIX}")
    message("Add benchmark: ${executable}")
    add_executable(${executable} matmul_ompss2.c)
    target_compile_definitions(${executable}
      PUBLIC -DISSTRONG=${STRONGTYPE}
      PUBLIC -DISMATVEC=${ISMATVEC})
    target_link_libraries(${executable} argparser)
  endforeach ()

  configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_${TEST}_dim.sh @ONLY)
  configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_${TEST}_dim.sh @ONLY)

endforeach ()

# install (copy) scripts
configure_file(validate.py validate.py @ONLY)

# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/create_trace.py" create_trace.py COPYONLY)

# Process
configure_file("${PROJECT_SOURCE_DIR}/utils/process_mn.py" process_mn.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/process_dim.py" process_dim.py COPYONLY)

