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


#SBATCH --qos=debug
#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

source argparse.sh
add_argument -a x -l exe -h "Executable file" -t file
add_argument -a d -l dim -h "Dimension" -t int
add_argument -a b -l bsize -h "Block size" -t int
add_argument -a i -l iter -h "Repetitions default[1]" -t int
add_argument -a o -l output -h "Output directory (used only for extrae)" -d ""
add_argument -a n -l nanos -h "Nanos version to use"
add_argument -a s -l scheduler -h "Nanos scheduler" -t string
parse_args "$@"
printargs >&2

env

module purge
module load bsc/1.0
module load gcc/7.2.0
module load liberep/git
module load mkl/2018.1
module load hwloc/1.11.8
module load impi/2018.1
module load boost/1.66.0
module load nanos6/cluster
module load mcxx/git

ulimit -c unlimited

# Environment variables
export NANOS6_CPU_SCHEDULER=fifo
export NANOS6_COMMUNICATION=mpi-2sided
export NANOS6_DISTRIBUTED_MEMORY=4G
export NANOS6_LOCAL_MEMORY=1G

export NANOS6_SCHEDULER=${ARGS[s]}
export NANOS6=${ARGS[n]}

iterations=${ARGS[i]}

# If EXTRAE
if [[ ${NANOS6} = extrae ]]; then

	module load EXTRAE/3.5.2
	dest="${ARGS[o]}/${SLURM_JOB_NAME}"
	export EXTRAE_CONFIG_FILE="${dest}.xml"
	iterations=1  # with extrae enabled only 1 iteration is allowed

	# Create xml for this for this benchmark if extrae is set
	sed -e "s|PREFIX|${SLURM_JOB_NAME}|"               \
		-e "s|EXTRAEHOME|${EXTRAE_HOME}|"              \
		-e "s|TMPDIR|${dest}_tmp|"\
		-e "s|FINALDIR|${ARGS[o]}|"  \
		-e "s|BINARY|${ARGS[x]}|"                      \
		-e "s|NPR|${SLURM_JOB_NUM_NODES}|"             \
		extrae_template.xml > ${EXTRAE_CONFIG_FILE}

	export NANOS6_EXTRAE_AS_THREADS=1
    export LD_PRELOAD="${EXTRAE_HOME}/lib/libnanosmpitrace.so"

	echo -e "# Enabled EXTRAE with EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE}"
	echo -e "# LD_PRELOAD=${LD_PRELOAD}"
fi

export NANOS6_COMMUNICATION=mpi-2sided
export NANOS6_DISTRIBUTED_MEMORY=4G
export NANOS6_LOCAL_MEMORY=1G

# Start run here
echo -e "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo -e "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo -e "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo -e "# QOS: ${SLURM_JOB_QOS}"
echo -e "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
echo -e "# Command: ${ARGS[x]} ${ARGS[d]} ${ARGS[b]} 0"
echo -e "# Iterations: ${iterations}"
env | grep NANOS6 | sed -e 's/^#*/# /'
echo -e "# --------------------------------------\n"

for ((it=0; it<iterations; ++it)) {
		echo "# Starting it: ${it} at: " $(date)
		start=$SECONDS
		srun taskset -c 0-3 ${ARGS[x]} ${ARGS[d]} ${ARGS[b]} 0
		end=$SECONDS
		echo "# Ending: " $(date)
		echo "# Elapsed: "$((end-start))
		echo -e "\n# --------------------------------------\n"
	}
