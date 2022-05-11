#!/bin/bash

# Copyright (C) 2022  Jimmy Aguilar Mena

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

source @PROJECT_BINARY_DIR@/argparse.sh

# Arguments
add_argument -a R -l repeats -h "Program repetitions default[3]" -t int -d 3
add_argument -a N -l ntasks -h "Number of tasks (mpi ranks)" -t list -d 1,2,4,8,16,32

add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program interations default[5]" -t int -d 5

parse_args "$@"

echo "# Execution start time: $(date)"
echo "# Args: $*"
printargs "# "

# Briefly check the executable validity as our argparser can't know.
if [ -z "${ARGS[REST]}" ]; then  # Check that there are executables to run
    echo "Error: No input executable provided ('ARGS[REST]' is empty)" >&2
    exit 1
fi

for EXE in ${ARGS[REST]}; do
    if [ ! -x $EXE ]; then       # Check all the files are executable
        echo "Error: '$EXE' is not an executable file" >&2
        exit 2
    fi

    PREFIX=${EXE%%_*}
    if [ -z ${TEST} ]; then      # check all inputs have a common prefix.
        TEST=${PREFIX}
    elif [ ${TEST} != ${PREFIX} ]; then
        echo "Error: No common prefix (${TEST} != ${PREFIX})" >&2
        exit 3
    fi
done

EXES=${ARGS[REST]}
EXES_COUNT=$(echo $EXES | wc -w)

# If we are here then we can start execution.
for NTASKS in ${ARGS[N]}; do
	for EXE in ${EXES}; do
		echo "# Starting executable: ${EXE} $((++EXECOUNT))/${EXES_COUNT}"

		COMMAND="mpirun -np=${NTASKS} ./${EXE} ${ARGS[D]} ${ARGS[B]} ${ARGS[I]}"

		echo -e "# Starting command: ${COMMAND}"
		echo "# ======================================"
		for it in $(seq ${REPEATS}); do
			echo "# Starting it: ${it} $(date)"
			start=${SECONDS}
			${COMMAND}
			end=${SECONDS}
			echo "# Ending: $(date)"
			echo "# Elapsed: $((end-start)) accumulated $((end-init))"
			echo "# --------------------------------------"
		done
	done
done # NTASKS
echo ""
