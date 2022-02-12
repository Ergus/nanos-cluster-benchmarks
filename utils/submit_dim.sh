#!/bin/bash

#SBATCH --signal=B:SIGUSR1@120

function got_signal() {
	echo "# $(date): Got exit signal"
}

trap got_signal SIGUSR1

# Declare command line arguments.
# In this script there are not default values to prevent errors.
source @PROJECT_BINARY_DIR@/argparse.sh

# Arguments to use in this script
add_argument -a R -l repeats -h "Repetitions per program" -t int
add_argument -a N -l ntasks -h "number of tasks" -t int
add_argument -a C -l cores -h "Number of cores per task" -t int

# Arguments for the executable
add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program interations" -t int

# Parse input command line arguments
parse_args "$@"

init=${SECONDS}

DIM=${ARGS[D]}
BS=${ARGS[B]}
REPEATS=${ARGS[R]}
ITS=${ARGS[I]}

NTASKS=${ARGS[N]}
CORES=${ARGS[C]}

# special nanos variables needed to set.
export NANOS6_CONFIG=@PROJECT_BINARY_DIR@/nanos6.toml

# Start run here printing run info header
echo "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo "# Nodes: ${SLURM_JOB_NUM_NODES} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo "# Ntasks: ${NTASKS} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE}"
echo "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo "# QOS: ${SLURM_JOB_QOS}"
echo "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
echo "# Walltime: $(squeue -h -j $SLURM_JOBID -o "%l")"

# Print command line arguments
printargs "# "

# Print nanos6 environment variables
env | grep NANOS6 | sed -e 's/^#*/# /'
echo ""

if [ $((SLURM_JOB_NUM_NODES*BS<=DIM)) != 1 ]; then
	echo "# Jump combination nodes: $node, dim: $rows, bs: $BS"
	exit
fi

EXES=${ARGS[REST]}
EXES_COUNT=$(echo $EXES | wc -w)

for EXE in ${EXES}; do
	# Check that EXE exists and is an executable
	if [[ ! -x $EXE ]]; then
		echo "# WARN: Input: ${EXE} is not an executable file"
		continue
	fi

	echo "# Starting executable: ${EXE} $((++EXECOUNT))/${EXES_COUNT}"
	for DISABLE_REMOTE in false true; do  # namespace enable/disable

		# Run only once for mpi benchmarks
		[ $DISABLE_REMOTE == true ] && [ ${EXE##*_} == 'mpi' ] && continue

		export NANOS6_CONFIG_OVERRIDE="cluster.disable_remote=${DISABLE_REMOTE}"

		COMMAND="srun --cpu-bind=cores --ntasks=${NTASKS} --cpus-per-task=${CORES} ./${EXE} $DIM $BS $ITS"

		echo -e "# Starting command: ${COMMAND}"
		echo "# ======================================"
		for it in $(seq ${REPEATS}); do
			echo "# Starting it: ${it} $(date)"
			start=${SECONDS}
			${COMMAND}
			end=${SECONDS}
			echo "# Ending: $(date)"
			echo "# Elapsed: $((end-start)) accumulated $((end-init))"
			echo "# --------------------------------------"
		done
	done
done

# We arrive here only when not wall time was reached.
# else the function job_finish_hook is called.
finalize=${SECONDS}
echo "# Done:  $(date)"
echo "# Elapsed: $((finalize-init))"
echo "# ======================================"
