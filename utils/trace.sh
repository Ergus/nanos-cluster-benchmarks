#!/usr/bin/env bash

export EXTRAE_ON=1
export EXTRAE_CONFIG_FILE=extrae.xml

export NANOS6_CONFIG=@PROJECT_BINARY_DIR@/nanos6.toml
export NANOS6_CONFIG_OVERRIDE="version.debug=false,version.instrument="extrae",cluster.disable_remote=false"

# Remove previous traces
if  [ -z ${MPI_LOCALRANKID} ] || [ $MPI_LOCALRANKID -eq 0 ]; then
	echo "Deleting old traces"
	rm -rf TRACE.* set-0/ *.prv *.row *.pcf
fi

export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so
$@
unset LD_PRELOAD

# Create the mpits
if  [ -z ${MPI_LOCALRANKID} ] || [ $MPI_LOCALRANKID -eq 0 ]; then
	echo "Creating mpits file manually"
	ls ${PWD}/set-0/*.mpit | sed 's/$/ named/' > TRACE.mpits
fi
