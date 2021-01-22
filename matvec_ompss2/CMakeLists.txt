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
add_executable(matvec_strong matvec_strong.cpp)
target_link_libraries (matvec_strong argparser)

add_executable(matvec_weak_nofetch matvec_weak.cpp)
target_link_libraries(matvec_weak_nofetch argparser)
target_compile_definitions(matvec_weak_nofetch PUBLIC -DFETCHTASK=0)

add_executable(matvec_weak_fetchall matvec_weak.cpp)
target_link_libraries(matvec_weak_fetchall argparser)
target_compile_definitions(matvec_weak_fetchall PUBLIC -DFETCHTASK=1)

add_executable(matvec_weak_fetchfirst matvec_weak.cpp)
target_link_libraries(matvec_weak_fetchfirst argparser)
target_compile_definitions(matvec_weak_fetchfirst PUBLIC -DFETCHTASK=2)

# install (copy) scripts
configure_file(validate.py validate.py @ONLY)
configure_file(submiter_mn.sh submiter_mn.sh @ONLY)
configure_file(submit_mn.sh submit_mn.sh @ONLY)
configure_file(extrae_template.xml extrae_template.xml COPYONLY)
configure_file(trace.sh trace.sh COPYONLY)
