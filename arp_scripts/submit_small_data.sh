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
		-F|--full
			If specified, translate everything
		-i|--image
			If specified, translate everything & save area detectors as images
		-t|--test
			Run the slurm job as test only to get job info
		-b|--bsub
			Run the job through LFS
		-E|--epicsAll
			If specified, translate all EPIVS PVs
EOF

}
#Look for args we expect, for now ignore other args
#Since we can have a mix of flags/args and args do in loop
for arg in "$@" 
do
	case $arg in
		-h|--help)
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
		-l|--locally)
			LOCALLY=true
			shift
			;;
		-f|--fast)
			NORECORDER=true
			shift
			;;
		-F|--full)
			FULL=True
			shift
			;;
		-i|--image)
			IMAGE=True
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
		-E|--epicsAll)
			EPICSALL=True
			shift
			;;
	esac
done

#Define cores if we don't have them
#Set to 1 if single is set
CORES=${CORES:='1'}
if [[ "$SINGLE" = true ]]; then
	CORES=1
fi
ARGS=''
#Generate smalldata_producer_arp.py args, for now if
#Experiment and Run not specified assume it's being handed
#by the ARP args in the os.environment
if [[ -v RUN ]]; then
	ARGS+=' --run '$RUN
fi
if [[ -v EXP ]]; then
	ARGS+=' --exp '$EXP
fi
if [[ -v NEVENTS ]]; then
	ARGS+=' --nevents '$NEVENTS
fi
if [[ -v DIRECTORY ]]; then
	ARGS+=' --dir '$DIRECTORY
fi
if [[ -v OFFLINE ]]; then
	ARGS+=' --offline'
fi
if [[ -v GATHER_INTERVAL ]]; then
	ARGS+=' --gather '$GATHER_INTERVAL
fi
if [[ -v NORECORDER ]]; then
	ARGS+=' --norecorder'
fi
if [[ -v EPICSALL ]]; then
	ARGS+=' --epicsAll'
fi
if [[ -v FULL ]]; then
	ARGS+=' --full'
fi
if [[ -v IMAGE ]]; then
	ARGS+=' --image'
fi

source /reg/g/psdm/etc/psconda.sh
ABS_PATH=/reg/g/psdm/sw/tools/smalldata_tools/examples
if [[ -v LOCALLY ]]; then
    #this should ideally be a relative path like ../smalldata_tools 
    #that did not work for me though....
    ABS_PATH=/reg/neh/home4/snelson/git_smd/smalldata_tools/examples
fi
if [[ -v SINGLE ]]; then
    $ABS_PATH/smalldata_producer_arp.py $ARGS
#else
#    sbatch --cpus-per-task=$CORES -p anagpu $ABS_PATH/smalldata_producer_arp.py $ARGS
fi
