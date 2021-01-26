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
add_argument -a x -l exe -h "Executable file" -t file
add_argument -a w -l wtime -h "Wall time limit for jobs" -t timer -d 01:00:00
add_argument -a q -l queue -h "queue" -t enum -e "debug bsc_cs xlarge" -d "bsc_cs"
add_argument -a R -l repeats -h "Repetitions per program default[1]" -t int -d 5

parse_args "$@"
printargs

now=$(date +%F_%H-%M-%S)
name=$(basename ${ARGS[x]})
resdir="results/${name}_${now}"

mkdir -p ${resdir}

echo "Output: ${resdir}"

nodes=(1 2 4 8)
echo "nodes: ${nodes[*]}"

for node in ${nodes[@]}; do
	echo "Submitting nodes: ${node}"

	jobname="${name}_${node}"
	filename="${resdir}/${jobname}"

 	sbatch --ntasks=${node} \
		   --time=${ARGS[w]} \
		   --qos=${ARGS[q]} \
 		   --job-name=${jobname} \
 		   --output="${resdir}/%x_%2a_%j.out" \
 		   --error="${resdir}/%x_%2a_%j.err" \
 		   ./submit_mn.sh -R ${ARGS[R]} -x ${ARGS[x]}
done
