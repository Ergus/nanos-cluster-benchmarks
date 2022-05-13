#!/bin/bash

# Copyright (C) 2021  Jimmy Aguilar Mena

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

# Arguments for this script
add_argument -a w -l wtime -h "Wall time limit for jobs" -t timer -d 12:00:00
add_argument -a q -l queue -h "Cluster queue" -t enum -e debug,bsc_cs,xlarge -d bsc_cs
add_argument -a o -l output -h "Output directory" -t string -d "results"

# Arguments to bypass
add_argument -a R -l repeats -h "Program repetitions default[3]" -t int -d 3
add_argument -a N -l ntasks -h "Number of tasks" -t list -d 1,2,4,8,16,32
add_argument -a C -l cores -h "Number of cores per node" -t list -d 12,24,48
add_argument -a F -l fitnodes -h "Fit processes in nodes" -t bool

add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t list
add_argument -a I -l iterations -h "Program interations default[5]" -t int -d 5

parse_args "$@"

mkdir -p "${ARGS[o]}"
logfile="${ARGS[o]}/submit.log"

# Redirect all the outputs to the logfile
exec &> >(tee -a ${ARGS[o]}/submit.log)

echo "# Submit time: $(date)"
echo "# Args: $*"
printargs "# "

# Check that there are executables to run
if [ -z "${ARGS[REST]}" ]; then
    echo "Error: No input executable provided ('ARGS[REST]' is empty)" >&2
    exit 1
fi

for EXE in ${ARGS[REST]}; do
    if [ ! -x $EXE ]; then    # Check all the files are executable
        echo "Error: '$EXE' is not an executable file" >&2
        exit 2
    fi

    PREFIX=${EXE%%_*}
    if [ -z ${TEST} ]; then   # check all inputs have a common prefix.
        TEST=${PREFIX}
    elif [ ${TEST} != ${PREFIX} ]; then
        echo "Error: No common prefix (${TEST} != ${PREFIX})" >&2
        exit 3
    fi
done

# If we are here then we can start the submission.
for CORES in ${ARGS[C]}; do
    for BS in ${ARGS[B]}; do

        JOBPREFIX="${TEST}_${ARGS[D]}_${BS}_${ARGS[I]}_${CORES}"

        OUTDIR="${ARGS[o]}/${JOBPREFIX}"
        mkdir -p ${OUTDIR}
        echo "# Output directory: ${OUTDIR}"

        for NTASKS in ${ARGS[N]}; do
            if [ ${ARGS[F]} == true ]; then # One process per node
                nodes=$(( (NTASKS * CORES + 48 - 1) / 48 ))
            else
                nodes=${NTASKS}
            fi
            command="sbatch --nodes=${nodes} \
                            --exclusive \
                            --time=${ARGS[w]} \
                            --qos=${ARGS[q]} \
                            --job-name="${JOBPREFIX}/${NTASKS}" \
                            --output="${ARGS[o]}/%x_%j.out" \
                            --error="${ARGS[o]}/%x_%j.err" \
                            --workdir=. \
                            ./submit_dim.sh \
                            -R ${ARGS[R]} \
                            -D ${ARGS[D]} \
                            -B ${BS} \
                            -I ${ARGS[I]} \
                            -N ${NTASKS}   \
                            -C ${CORES} \
                            ${ARGS[REST]} "

            # Print and execute command
            echo ${command// +/ }
            ${command}
        done # NTASKS
    done # BS
done # CORES
echo ""
