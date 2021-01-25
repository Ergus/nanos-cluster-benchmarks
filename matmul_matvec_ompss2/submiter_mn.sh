#!/bin/bash

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

source @PROJECT_BINARY_DIR@/argparse.sh
add_argument -a x -l exe -h "Executable file" -t file
add_argument -a i -l iter -h "Inner Repetitions default[1]" -t int -d 1
add_argument -a p -l print -h "Print matrices default[true]" -t bool -d 1
add_argument -a I -l Iter -h "Repetitions per program default[1]" -t int -d 1
add_argument -a w -l wtime -h "Wall time limit for jobs" -t timer -d 00:10:00
add_argument -a q -l queue -h "queue" -t enum -e "debug bsc_cs xlarge" -d "debug"

add_argument -a N -l Nodes -h "List of nodes numbers" -t string -d 4
add_argument -a d -l dimensions -h "List of dimensions" -t string -d 128
add_argument -a b -l bsizes -h "List block sizes" -t string -d 8

add_argument -a t -l type -h "Type of scalability experiment to test" -t enum\
	     -d "strong" -e "weak strong"

add_argument -a n -l nanos -h "Nanos version to use" -t enum\
	     -e "Debug optimized extrae verbose-debug" -d "debug"
add_argument -a s -l scheduler -h "Nanos scheduler" -t enum\
	     -e "cluster-random cluster-locality" -d "cluster-random"
add_argument -a m -l mpi -h "MPI version to use" -t enum\
	     -e "impi ompi" -d "impi"

#-e "distributed distributedLocation distributedLocationBalanced"\
parse_args "$@"
printargs

IFS=', ' read -r -a nodes <<< "${ARGS[N]}"
IFS=', ' read -r -a blocksize <<< "${ARGS[b]}"
IFS=', ' read -r -a dims <<< "${ARGS[d]}"

now=$(date +%F_%H-%M-%S)
name=${ARGS[x]/%.nanos6}
resdir="results/${name}_${ARGS[n]}_${ARGS[s]}_${ARGS[t]}_${ARGS[p]}_${ARGS[i]}_${now}"

mkdir -p ${resdir}

echo "Directory: ${resdir}"
echo "Submitting: group"
echo "nodes: ${nodes[*]}"
echo "blocks: ${blocksize[*]}"
echo "dim: ${dims[*]}"

for dim in ${dims[@]}; do
	for node in ${nodes[@]}; do

		case ${ARGS[t]} in
			strong) rows=${dim};;
			weak) rows=$(( dim * node ));;
			*) echo "Error in experiment type";;
		esac

		for bs in ${blocksize[@]}; do

			if [ $((node*bs<=dim)) = 1 ]; then
				printf "Submitting combination dim: %s bs: %s nodes: %s\n" \
					   $dim $bs $node
				jobname="mv_${dim}_${bs}_${node}"
				filename="${resdir}/${jobname}"

 				sbatch --ntasks=${node} \
				       --time=${ARGS[w]} \
				       --qos=${ARGS[q]} \
 				       --job-name=${jobname} \
 				       --output="${resdir}/%x_%2a_%j.out" \
 				       --error="${resdir}/%x_%2a_%j.err" \
 				       ./submit_mn.sh -x ${ARGS[x]} \
				       -R ${rows} -C ${dim} -b ${bs} \
				       -i ${ARGS[i]} \
				       -p $([[ ${ARGS[p]} = true ]] && echo 1 || echo 0) \
				       -I ${ARGS[I]} -o ${resdir} \
				       -n ${ARGS[n]} -s ${ARGS[s]} -m ${ARGS[m]}
			else
				printf "Jump combination dim: %s bs: %s nodes: %s\n" \
				       $dim $bs $node
			fi
		done
	done
done
