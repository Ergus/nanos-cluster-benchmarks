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
add_argument -a N -l ntasks -h "number of tasks" -t list
add_argument -a C -l cores -h "Number of cores per task" -t int

# Arguments for the executable
add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program interations" -t int

add_argument -a M -l disable_remote -h "Disable reMote" -t enum -e true,false
add_argument -a L -l leader -h "Leader thread" -t enum -e true,false
add_argument -a G -l group -h "Group Messages" -t enum -e true,false
add_argument -a W -l writeid -h "Use Writeid" -t enum -e true,false
add_argument -a H -l helpers -h "Number of Handler helpers" -t int

# Parse input command line arguments
parse_args "$@"

init=${SECONDS}

DIM=${ARGS[D]}
BS=${ARGS[B]}

# special nanos variables needed to set.
export NANOS6_CONFIG=@PROJECT_BINARY_DIR@/nanos6.toml

NANOS6_CONFIG_OVERRIDE="cluster.disable_remote=${ARGS[M]}"
NANOS6_CONFIG_OVERRIDE+=",cpumanager.reserve_cpu_leader=${ARGS[L]}"
NANOS6_CONFIG_OVERRIDE+=",cluster.group_messages=${ARGS[G]}"
NANOS6_CONFIG_OVERRIDE+=",cluster.enable_write_id=${ARGS[W]}"
NANOS6_CONFIG_OVERRIDE+=",cluster.num_message_handler_workers=${ARGS[H]}"

export NANOS6_CONFIG_OVERRIDE=${NANOS6_CONFIG_OVERRIDE}

# Start run here printing run info header
echo "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo "# Nodes: ${SLURM_JOB_NUM_NODES} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo "# Ntasks: ${ARGS[N]} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE}"
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

	for NTASKS in ${ARGS[N]}; do
		COMMAND="srun --cpu-bind=cores --ntasks=${NTASKS} --cpus-per-task=${ARGS[C]} ./${EXE} $DIM $BS ${ARGS[I]}"

		echo -e "# Starting command: ${COMMAND}"
		echo "# ======================================"
		for it in $(seq ${ARGS[R]}); do
			echo "# Starting it: ${it} $(date)"
			start=${SECONDS}
			${COMMAND}
			end=${SECONDS}
			echo "# Ending: $(date)"
			echo "# Elapsed: $((end-start)) accumulated $((end-init))"
			echo "# --------------------------------------"
		done
	done # ntasks
done

# We arrive here only when not wall time was reached.
# else the function job_finish_hook is called.
finalize=${SECONDS}
echo "# Done:  $(date)"
echo "# Elapsed: $((finalize-init))"
echo "# ======================================"
