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

add_executable(matvec_weak_nofetch matvec.c)
target_compile_definitions(matvec_weak_nofetch
  PUBLIC -DFETCHTASK=0
  PUBLIC -DISMATVEC=1)
target_link_libraries(matvec_weak_nofetch argparser)

add_executable(matvec_weak_fetchall matvec.c)
target_compile_definitions(matvec_weak_fetchall
  PUBLIC -DFETCHTASK=1
  PUBLIC -DISMATVEC=1)
target_link_libraries(matvec_weak_fetchall argparser)

add_executable(matvec_weak_fetchfirst matvec.c)
target_compile_definitions(matvec_weak_fetchfirst
  PUBLIC -DFETCHTASK=2
  PUBLIC -DISMATVEC=1)
target_link_libraries(matvec_weak_fetchfirst argparser)

# Add the executables and link
add_executable(matmul_strong matvec.c)
target_compile_definitions(matmul_strong
  PUBLIC -DISSTRONG
  PUBLIC -DISMATVEC=0)
target_link_libraries (matmul_strong argparser)

add_executable(matmul_weak_nofetch matvec.c)
target_compile_definitions(matmul_weak_nofetch
  PUBLIC -DFETCHTASK=0
  PUBLIC -DISMATVEC=0)
target_link_libraries(matmul_weak_nofetch argparser)

add_executable(matmul_weak_fetchall matvec.c)
target_compile_definitions(matmul_weak_fetchall
  PUBLIC -DFETCHTASK=1
  PUBLIC -DISMATVEC=0)
target_link_libraries(matmul_weak_fetchall argparser)

add_executable(matmul_weak_fetchfirst matvec.c)
target_compile_definitions(matmul_weak_fetchfirst
  PUBLIC -DFETCHTASK=2
  PUBLIC -DISMATVEC=0)
target_link_libraries(matmul_weak_fetchfirst argparser)


# install (copy) scripts
configure_file(validate.py validate.py @ONLY)
configure_file(submiter_mn.sh submiter_mn.sh @ONLY)
configure_file(submit_mn.sh submit_mn.sh @ONLY)
configure_file(extrae_template.xml extrae_template.xml COPYONLY)
configure_file(trace.sh trace.sh COPYONLY)