#!/bin/bash

#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --signal=B:USR1@10

# Some general utils for when killed by time.
init=${SECONDS}
job_finish_hook()
{
	finalize=${SECONDS}
	echo "# Killing: $(date)"
	echo "# Elapsed: $((finalize-init))"
	echo "# ======================================"
}
trap 'job_finish_hook' USR1

# Declare command line arguments.
source @PROJECT_BINARY_DIR@/argparse.sh
add_argument -a x -l exe -h "Executable file" -t file
add_argument -a R -l repeats -h "Repetitions per program" -t int
add_argument -a N -l namespace -h "Namespace propagation enabled" -t int
add_argument -a I -l iterations -h "Program interations" -t int

# Parse input command line arguments
parse_args "$@"

dims=(16384 32768 65536)
blocksizes=(8 16 32 64 128 256 512)

REPEATS=${ARGS[R]}
ITS=${ARGS[I]}

# special nanos variables needed to set.
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

# Print command line arguments
printargs "# "

# Print nanos6 environment variables
env | grep NANOS6 | sed -e 's/^#*/# /'
echo "# ======================================"

for dim in ${dims[@]}; do
	for bs in ${blocksizes[@]}; do
		if [ $((SLURM_JOB_NUM_NODES*bs<=dim)) = 1 ]; then
			COMMAND="${ARGS[x]} $dim $bs $ITS"
			echo -e "# Starting command: ${COMMAND}"
			echo "# ======================================"
			for ((it=0; it<${REPEATS}; ++it)) {
				echo "# Starting it: ${it} at: $(date)"
				start=${SECONDS}
				srun ${COMMAND}
				end=${SECONDS}
				echo "# Ending:  $(date)"
				echo "# Elapsed: $((end-start))"
				echo "# --------------------------------------"
			}
			echo ""
		else
			echo "# Jump combination nodes: $node, dim: $dim, bs: $bs"
		fi
	done
done

# We arrive here only when not wall time was reached.
# else the function job_finish_hook is called.
finalize=${SECONDS}
echo "# Done:  $(date)"
echo "# Elapsed: $((finalize-init))"
echo "# ======================================"
