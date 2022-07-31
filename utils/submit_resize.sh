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
add_argument -a R -l repeats -h "Repetitions" -t int
add_argument -a N -l ntasks -h "Number of tasks (sizes)" -t list
add_argument -a C -l cores -h "Number of cores per task" -t int

# Arguments for the executable
add_argument -a D -l dim -h "Matrix dimension" -t int
add_argument -a B -l BS -h "Blocksize" -t int
add_argument -a I -l iterations -h "Program iterations" -t int

# Parse input command line arguments
parse_args "$@"

init=${SECONDS}

DIM=${ARGS[D]}
BS=${ARGS[B]}
REPEATS=${ARGS[R]}
ITS=${ARGS[I]}

NTASKS=${ARGS[N]}
CORES=${ARGS[C]}

# The REST argument contains the list of executables...

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

# List of nodes growing and decreasing
GROW=$(echo -e ${ARGS[N]// /\\n} | sort -n | tr '\n' ' ')
FALL=$(echo -e ${ARGS[N]// /\\n} | sort -n -r | tr '\n' ' ')

MINNTASK=$(echo -e ${GROW// /\\n} | head -n1)
MAXNTASK=$(echo -e ${FALL// /\\n} | head -n1)

export NANOS6_CONFIG_OVERRIDE="cluster.num_max_nodes=${MAXNTASK}"

for EXE in ${EXES}; do
	# Check that EXE exists and is an executable
	if [[ ! -x $EXE ]]; then
		echo "# WARN: Input: ${EXE} is not an executable file"
		continue
	fi

	echo "# Starting executable: ${EXE} $((++EXECOUNT))/${EXES_COUNT}"
	echo "# ======================================"
	if [[ ${EXE} =~ "_ompss2" ]]; then
		for IT in $(seq ${REPEATS}); do
			for NTASKS in ${ARGS[N]}; do
				echo "# Starting it: ${it} $(date)"
				start=${SECONDS}

				srun --cpu-bind=cores --ntasks=${NTASKS} --cpus-per-task=${CORES} \
					 ./${EXE} $DIM $BS $ITS 0 $((IT + 1)) ${IT} $(date "+%s %N")

				[[ ${IT} == 0 ]] || rm -r ${IT}

				end=${SECONDS}
				echo "# Ending: $(date)"
				echo "# Elapsed: $((end-start)) accumulated $((end-init))"
				echo "# --------------------------------------"
			done
		done
	elif [[ ${EXE} =~ "_resize" ]]; then
		echo "# Starting it: $(date)"
		start=${SECONDS}

		srun --cpu-bind=cores --ntasks=${MINNTASK} --cpus-per-task=${CORES} \
			 ./${EXE} ${DIM} ${BS} ${GROW} ${FALL}

		end=${SECONDS}
		echo "# Ending: $(date)"
		echo "# Elapsed: $((end-start)) accumulated $((end-init))"
		echo "# --------------------------------------"
	fi

done

# We arrive here only when not wall time was reached.
# else the function job_finish_hook is called.
finalize=${SECONDS}
echo "# Done:  $(date)"
echo "# Elapsed: $((finalize-init))"
echo "# ======================================"
