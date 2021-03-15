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

SET(FETCHLIST nofetch fetchall fetchfirst)
foreach (FETCHTASK 0;1;2)
  list(GET FETCHLIST ${FETCHTASK} FETCHSUFIX)

  SET(NODELIST nonode node)

  set(executable "${TEST}_${FETCHSUFIX}")
  message("Add benchmark: ${executable}")
  add_executable(${executable} jacobi.c)
  target_compile_definitions(${executable}
    PUBLIC -DFETCHTASK=${FETCHTASK})
  target_link_libraries(${executable} argparser)

endforeach ()

# install (copy) scripts
# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_${TEST}_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_${TEST}_dim.sh @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

#Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/create_trace.py" create_trace.py COPYONLY)

# Process
configure_file("${PROJECT_SOURCE_DIR}/utils/process_mn.py" process_mn.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/process_dim.py" process_dim.py COPYONLY)

