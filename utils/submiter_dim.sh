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
add_argument -a C -l cores -h "Number of cores per node" -t int -d 48

add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program interations default[5]" -t int -d 5

parse_args "$@"
printargs "# "

jobname="@TEST@_${ARGS[D]}_${ARGS[B]}_${ARGS[I]}_${ARGS[C]}"

mkdir -p "results/${jobname}"
echo "# Output directory: ${jobname}"

ntasks=(1 2 4 8 16 32)
echo "# List num ntasks: ${nodes[*]}"

for ntask in ${ntasks[@]}; do
	echo "# Submitting for ${ntask} task[s]"

 	sbatch --ntasks=${ntask} \
		   --time=${ARGS[w]} \
		   --qos=${ARGS[q]} \
 		   --job-name="${jobname}/${ntask}" \
 		   --output="results/%x_%j.out" \
 		   --error="results/%x_%j.err" \
 		   ./submit_@TEST@_dim.sh \
		   -R ${ARGS[R]} \
		   -D ${ARGS[D]} \
		   -B ${ARGS[B]} \
		   -I ${ARGS[I]} \
		   -N ${ntask}   \
		   -C ${ARGS[C]}
done
