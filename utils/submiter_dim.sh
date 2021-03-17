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
add_argument -a q -l queue -h "Cluster queue" -t enum -e "debug bsc_cs xlarge" -d "bsc_cs"

# Arguments to bypass
add_argument -a R -l repeats -h "Program repetitions default[3]" -t int -d 3
add_argument -a W -l weakscaling -h "Do weak scaling (dim*sqrt(nodes))" -t int -d 0

add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program interations default[5]" -t int -d 5

parse_args "$@"
printargs "# "

[ ${ARGS[W]} = 1 ] && suffix="weak" || suffix="strong"
resdir="results/@TEST@_${ARGS[D]}_${ARGS[B]}_${ARGS[I]}_${suffix}"

mkdir -p ${resdir}
echo "# Output directory: ${resdir}"

nodes=(1 2 4 8 16)
echo "# List num nodes: ${nodes[*]}"

for node in ${nodes[@]}; do
	echo "# Submitting for ${node} node[s]"

	jobname="@TEST@_${node}"

 	sbatch --ntasks=${node} \
		   --time=${ARGS[w]} \
		   --qos=${ARGS[q]} \
 		   --job-name=${jobname} \
 		   --output="${resdir}/%x_%j.out" \
 		   --error="${resdir}/%x_%j.err" \
 		   ./submit_@TEST@_dim.sh -R ${ARGS[R]} -D ${ARGS[D]} -B ${ARGS[B]} -I ${ARGS[I]} -W ${ARGS[W]}
done
