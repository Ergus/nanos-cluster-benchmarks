#!/usr/bin/env bash

echo "Running with Extrae rank ${PMPI_RANK}"

export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so

taskset -c 0-2 $@
