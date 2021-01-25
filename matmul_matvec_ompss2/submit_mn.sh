#!/bin/bash

#SBATCH --workdir=.

#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48

source @PROJECT_BINARY_DIR@/argparse.sh
add_argument -a x -l exe -h "Executable file" -t file
add_argument -a R -l rows -h "Rows" -t int
add_argument -a C -l cols -h "Columns" -t int
add_argument -a b -l bsize -h "Block size" -t int
add_argument -a i -l iter -h "Inner Repetitions in task default[1]" -t int
add_argument -a p -l print -h "Print matrix default" -t int
add_argument -a I -l Iter -h "Repetitions per program default[1]" -t int
add_argument -a o -l output -h "Output directory (used only for extrae)" -d ""
add_argument -a n -l nanos -h "Nanos version to use"
add_argument -a s -l scheduler -h "Nanos scheduler" -t string
add_argument -a m -l mpi -h "MPI version" -t string
parse_args "$@"
printargs >&2

module purge
module load bsc/1.0
module load gcc/7.2.0

if [[ ${ARGS[m]} == "impi" ]]; then
	module load hwloc/1.11.8
	module load intel/2018.1 impi/2018.1 mkl/2018.1
	module load EXTRAE/3.7.1
	module load boost/1.64.0-mpi
	module load nanos6/cluster
	module load mcxx/git_impi
elif [[ ${ARGS[m]} == "ompi" ]]; then
	module load openmpi/3.1.1
	module load EXTRAE/3.6.1
	module load nanos6/cluster_ompi
	module load mcxx/git_ompi
else
	echo "ERROR: You need to specify impi or ompi as first argument"
fi


ulimit -c unlimited

# Environment variables
export NANOS6_CPU_SCHEDULER=fifo
export NANOS6_COMMUNICATION=mpi-2sided
export NANOS6_DISTRIBUTED_MEMORY=128G
export NANOS6_LOCAL_MEMORY=8G

export NANOS6_SCHEDULER=${ARGS[s]}
export NANOS6=${ARGS[n]}

iterations=${ARGS[I]}

COMMAND="${ARGS[x]} ${ARGS[R]} ${ARGS[C]} ${ARGS[b]} ${ARGS[i]} ${ARGS[p]}"

# If EXTRAE
if [[ ${NANOS6} = extrae ]]; then

	dest="${ARGS[o]}/${SLURM_JOB_NAME}"
	export EXTRAE_CONFIG_FILE="${dest}.xml"
	iterations=1  # with extrae enabled only 1 iteration is allowed

	# Create xml for this for this benchmark if extrae is set
	sed -e "s|PREFIX|${SLURM_JOB_NAME}|"    \
	    -e "s|EXTRAEHOME|${EXTRAE_HOME}|"   \
	    -e "s|TMPDIR|${dest}_tmp|"          \
	    -e "s|FINALDIR|${ARGS[o]}|"         \
	    -e "s|BINARY|${ARGS[x]}|"           \
	    -e "s|NPR|${SLURM_JOB_NUM_NODES}|"  \
	    extrae_template.xml > ${EXTRAE_CONFIG_FILE}

	export NANOS6_EXTRAE_AS_THREADS=1
	#export LD_PRELOAD="${EXTRAE_HOME}/lib/libnanosmpitrace.so"

	#echo -e "# Enabled EXTRAE with EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE}"
	#echo -e "# LD_PRELOAD=${LD_PRELOAD}"

	COMMAND="./trace.sh ${COMMAND}"
fi

# Start run here
echo -e "# Job: ${SLURM_JOB_NAME} id: ${SLURM_JOB_ID}"
echo -e "# Nodes: ${SLURM_JOB_NUM_NODES} Tasks_per_Node: ${SLURM_NTASKS_PER_NODE} Cores_per_node: ${SLURM_JOB_CPUS_PER_NODE}"
echo -e "# Nodes_List: ${SLURM_JOB_NODELIST}"
echo -e "# QOS: ${SLURM_JOB_QOS}"
echo -e "# Account: ${SLURM_JOB_ACCOUNT} Submitter_host: ${SLURM_SUBMIT_HOST} Running_Host: ${SLURMD_NODENAME}"
echo -e "# Iterations: ${iterations}"
echo -e "# Command: " ${COMMAND}
env | grep NANOS6 | sed -e 's/^#*/# /'
echo -e "# ======================================\n"

for ((it=0; it<iterations; ++it))
    {
	    echo "# Starting it: ${it} at: " $(date)
	    start=$SECONDS
	    srun ${COMMAND}
	    end=$SECONDS
	    echo "# Ending: " $(date)
	    echo "# Elapsed: "$((end-start))
	    echo -e "\n# --------------------------------------\n"
    }
