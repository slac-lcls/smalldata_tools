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
		-x|--norecorder
			If specified, don't use recorder data
		-F|--full
			eIf specified, translate everything
		-t|--test
			Run the slurm job as test only to get job info
EOF

}
#Look for args we expect, for now ignore other args
#Since we can have a mix of flags/args and args do in loop

#Generate smalldata_producer_arp.py args, for now if
#Experiment and Run not specified assume it's being handed
#by the ARP args in the os.environment

POSITIONAL=()
while [[ $# -gt 0 ]]
do
        key="$1"

	case $key in
		-h|--help)
			usage
			exit
			;;
		-l|--locally)
			LOCALLY=true
			shift
			;;
		-F|--full)
			FULL=True
			shift
			;;
                 *)
                        POSITIONAL+=("$1")
			shift
			;;                     
	esac
done
set -- "${POSITIONAL[@]}"

source /reg/g/psdm/etc/psconda.sh
ABS_PATH=/reg/g/psdm/sw/tools/smalldata_tools/examples
if [[ -v LOCALLY ]]; then
    ABS_PATH=/reg/d/psdm/xcs/xcsx43118/results/smalldata_tools/examples
fi

if [[ -v FULL ]]; then
    mpirun $ABS_PATH/smalldata_producer_full_arp.py $@
else
    mpirun $ABS_PATH/smalldata_producer_arp.py $@
fi
