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
add_argument -a N -l nodes -h "Number of nodes" -t list -d 1,2,4,8,16,32
add_argument -a C -l cores -h "Number of cores per node" -t list -d 12,24,48


add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t list
add_argument -a I -l iterations -h "Program interations default[5]" -t int -d 5

parse_args "$@"

mkdir -p "${ARGS[o]}"
logfile="${ARGS[o]}/submit.log"
echo "# Submit time: $(date)" | tee -a ${logfile}

printargs "# " | tee -a ${logfile}

if [ -z "${ARGS[REST]}" ]; then
    echo "Error: No input executable provided ('ARGS[REST]' is empty)" >&2
    exit 1
fi

for EXE in ${ARGS[REST]}; do
    if [ ! -x $EXE ]; then
        echo "Error: '$EXE' is not an executable file" >&2
        exit 2
    fi
done

echo "# List ntasks: [ ${ARGS[N]} ]"

for BS in ${ARGS[B]}; do

    jobname="@TEST@_${ARGS[D]}_${BS}_${ARGS[I]}"

    OUTDIR="${ARGS[o]}/${jobname}"
    mkdir -p ${OUTDIR}
    echo "# Output directory: ${OUTDIR}"

    for ntask in ${ARGS[N]}; do

        command="sbatch --ntasks=${ntask} \
                        --time=${ARGS[w]} \
                        --qos=${ARGS[q]} \
                        --job-name="${jobname}/${ntask}" \
                        --output="${ARGS[o]}/%x_%j.out" \
                        --error="${ARGS[o]}/%x_%j.err" \
                        ./submit_matvec_dim.sh \
                        -R ${ARGS[R]} \
                        -D ${ARGS[D]} \
                        -B ${BS} \
                        -I ${ARGS[I]} \
                        -N ${ntask}   \
                        -C ${ARGS[C]// /,} \
                        ${ARGS[REST]} "

        echo ${command// +/ }
        ${command}
    done # ntasks

done | tee -a ${logfile}
echo "" >> ${logfile}
