#!/usr/bin/env bash

export EXTRAE_ON=1
export EXTRAE_CONFIG_FILE=extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so

export NANOS6_CONFIG=@PROJECT_BINARY_DIR@/nanos6.toml
export NANOS6_CONFIG_OVERRIDE="version.debug=false,version.instrument="extrae",cluster.disable_remote=false"

./$@
