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
add_argument -a w -l wtime -h "Wall time limit for jobs" -t timer -d 06:00:00
add_argument -a q -l queue -h "Cluster queue" -t enum -e "debug bsc_cs xlarge" -d "bsc_cs"
add_argument -a R -l repeats -h "Program repetitions default[5]" -t int -d 5
add_argument -a N -l namespace -h "Namespace propagation enabled default[1]" -t int -d 1
add_argument -a I -l iterations -h "Program interations default[5]" -t int -d 5

parse_args "$@"
printargs "# "

name=$(basename ${ARGS[x]})
resdir="results/${name}_${ARGS[N]}"

mkdir -p ${resdir}
echo "# Output directory: ${resdir}"

nodes=(1 2 4 8 16)
echo "# List num nodes: ${nodes[*]}"

for node in ${nodes[@]}; do
	echo "# Submitting for ${node} node[s]"

	jobname="${name}_${node}"

 	sbatch --ntasks=${node} \
		   --time=${ARGS[w]} \
		   --qos=${ARGS[q]} \
 		   --job-name=${jobname} \
 		   --output="${resdir}/%x_%j.out" \
 		   --error="${resdir}/%x_%j.err" \
 		   ./submit_mn.sh -R ${ARGS[R]} -x ${ARGS[x]} -N ${ARGS[N]} -I ${ARGS[I]}
done
