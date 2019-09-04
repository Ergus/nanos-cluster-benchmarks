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

# This benchmark uses Liberep
find_package (LIBEREP REQUIRED)
include_directories (${LIBEREP_INCLUDE})

# Add the executables and link
add_executable (matvec_strong matvec_strong.cpp)
target_link_libraries (matvec_strong ${LIBEREP_LIBRARY})

# install (copy) scripts
configure_file(submiter_mn.sh submiter_mn.sh COPYONLY)
configure_file(submit_mn.sh submit_mn.sh COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/utils/argparse.sh argparse.sh COPYONLY)
