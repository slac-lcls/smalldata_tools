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
		-n|--nevents
			Number of events to analyze
		-l|--locally
			If specified, will run locally
		-x|--norecorder
			If specified, don't use recorder data
		-D|--default
			If specified, translate only small data
		-f|--full
			eIf specified, translate everything
                -i|--image
			If specified, translate everything & save area detectors as images
                -T|--tiff
			If specified, translate everything & save area detectors as images * single-event tiffs
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
		-e|--experiment)
                        POSITIONAL+=("--experiment $2")
			shift
			shift
			;;
		-r|--run)
                        POSITIONAL+=("--run $2")
			shift
			shift
			;;
		-n|--nevents)
		        NEVENTS=$2
                        POSITIONAL+=("--nevents $2")
			shift
			shift
			;;
		-d|--directory)
                        POSITIONAL+=("--directory $2")
			shift
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
    #DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
    ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/examples/g`
fi

if [ -v NEVENTS ] && [ $NEVENTS -lt 100 ]; then
    $ABS_PATH/smalldata_producer.py $@
else
    mpirun $ABS_PATH/smalldata_producer.py $@
fi
