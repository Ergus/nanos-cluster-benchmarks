#!/bin/bash

#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

# Example to use this script:
# sbatch --ntasks=2 --time=00:15:00 --qos=bsc_cs --job-name=ext_ch_2 ./submit_extrae.sh ./cholesky_weak 512 32

# Declare command line arguments.
source @PROJECT_BINARY_DIR@/argparse.sh
add_argument -a n -l np -h "Number of MPI processes" -t int -d ${SLURM_JOB_NUM_NODES}
add_argument -a N -l namespace -h "Namespace propagation enabled" -t int
parse_args "$@"

# Special nanos variables needed to set.
export NANOS6_CONFIG=@PROJECT_BINARY_DIR@/nanos6.toml

# Enable/disable namepase (parameter -N)
if [ ${ARGS[N]} = 1 ]; then
	export NANOS6_CONFIG_OVERRIDE="cluster.disable_remote=false"
else
	export NANOS6_CONFIG_OVERRIDE="cluster.disable_remote=true"
fi

# Start run here printing run info header
echo "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo "# QOS: ${SLURM_JOB_QOS}"
echo "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
echo "# Walltime: $(squeue -h -j $SLURM_JOBID -o "%l")"

# Print command line arguments
printargs "# "

# Print nanos6 environment variables
env | grep NANOS6 | sed -e 's/^#*/# /'
echo "# ======================================"

mpirun -np ${ARGS[n]} ./trace.sh ${ARGS[REST]}
