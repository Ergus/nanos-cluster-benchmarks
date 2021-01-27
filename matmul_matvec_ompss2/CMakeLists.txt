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
add_executable(matvec_strong matvec.c)
target_compile_definitions(matvec_strong
  PUBLIC -DISSTRONG
  PUBLIC -DISMATVEC=1)
target_link_libraries (matvec_strong argparser)

add_executable(matmul_strong matvec.c)
target_compile_definitions(matmul_strong
  PUBLIC -DISSTRONG
  PUBLIC -DISMATVEC=0)
target_link_libraries (matmul_strong argparser)


# The order is important
SET(TESTLIST matmul matvec)
foreach (ISMATVEC 0;1)
  list(GET TESTLIST ${ISMATVEC} TEST)

  SET(FETCHLIST nofetch fetchall fetchfirst)
  foreach (FETCHTASK 0;1;2)
    list(GET FETCHLIST ${FETCHTASK} FETCHSUFIX)

    SET(NODELIST nonode node)
    foreach (WITH_NODE 0;1)
      list(GET NODELIST ${WITH_NODE} NODESUFIX)

      set(executable "${TEST}_weak_${FETCHSUFIX}_${NODESUFIX}")
      message("Add benchmark: ${executable}")
      add_executable(${executable} matvec.c)
      target_compile_definitions(${executable}
	PUBLIC -DFETCHTASK=${FETCHTASK}
	PUBLIC -DWITH_NODE=${WITH_NODE}
	PUBLIC -DISMATVEC=${ISMATVEC})
      target_link_libraries(${executable} argparser)

    endforeach ()
  endforeach ()
endforeach ()

# install (copy) scripts
configure_file(submit_mn.sh submit_mn.sh @ONLY)
configure_file(submiter_mn.sh submiter_mn.sh @ONLY)

# configure_file(validate.py validate.py @ONLY)
# configure_file(extrae_template.xml extrae_template.xml COPYONLY)
# configure_file(trace.sh trace.sh COPYONLY)
