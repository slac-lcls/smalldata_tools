#!/bin/bash

# Took this from Silke's script, nice implementation
usage()
{
cat << EOF
$(basename "$0"): 
	Script to launch a smalldata_tools run analysis
	
	OPTIONS:
		-h|--help
			Definition of options
		-e|--experiment
			Experiment name (i.e. cxilr6716)
		-r|--run
			Run Number
		-d|--directory
			Full path to directory for output file
		-q|--queue
			Queue to use on SLURM
		-c|--cores
			Number of cores to be utilized
		-s|--single
			Run on a single core
		-n|--nevents
			Number of events to analyze
		-l|--locally
			If specified, will run locally
		-f|--fast
			If specified, don't use recorder data
		-t|--test
			Run the slurm job as test only to get job info
		-b|--bsub
			Run the job through LFS
EOF
}

#Look for args we expect, for now ignore other args
#Since we can have a mix of flags/args and args do in loop
for arg in "$@" 
do
	case $arg in
		-h| --help)
			usage
			exit
			;;
		-r|--run) 
			RUN="$2"
			shift
			shift
			;;
		-e|--experiment)
			EXP="$2"
			shift
			shift
			;;
		-d|--directory)
			DIRECTORY="$2"
			shift
			shift
			;;
		-q|--queue)
			QUEUE="$2"
			shift
			shift
			;;
		-c|--cores)
			CORES="$2"
			shift
			shift
			;;
		-s|--single)
			SINGLE=true
			shift
			;;
		-n|--nevents)
			NEVENTS="$2"
			shift
			shift
			;;
		-l|--local)
			LOCALLY=true
			shift
			;;
		-f|--fast)
			NORECORDER=true
			shift
			;;
		-o|--offline)
			OFFLINE=true
			shift
			;;
		-t|--test)
			TEST=true
			shift
			;;
		-b|--bsub)
			LFS=true
			shift
			;;
		-g|--gather_interval)
			GATHER_INTERVAL="$2"
			shift
			shift
			;;
	esac
done

#Activate Conda Environment
source /reg/g/psdm/etc/psconda.sh

#Path to script
ABS_PATH=/reg/g/psdm/sw/tools/smalldata_tools/examples

#Define cores if we don't have them
#Set to 1 if single is set
CORES=${CORES:='1'}
if [ "$SINGLE" = true ]; then
	CORES=1
fi

#Generate smalldata_producer_arp.py args
ARGS="--run ${RUN} "
ARGS+="--exp ${EXPERIMENT} "
if [[ -v NEVENTS ]]; then
	ARGS+="--nevt ${NEVENTS} "
fi
if [[ -v DIRECTORY ]]; then
	ARGS+="--dir ${DIRECTORY} "
fi
if [[ -v OFFLINE ]]; then
	ARGS+="--offline ${OFFLINE} "
fi
if [[ -v GATHER_INTERVAL ]]; then
	ARGS+="--gather ${GATHER_INTERVAL} "
fi
if [[ -v NORECORDER ]]; then
	ARGS+="--norecorder ${NORECORDER}"
fi

#sbatch --cpus-per-task=$CORES --test-only $ABS_PATH/smalldata_producer_arp.py

sbatch --nodes=1 --time=5 $ABS_PATH/smalldata_producer_arp.py "$ARGS"
