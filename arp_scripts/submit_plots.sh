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
		--experiment
			Experiment name (i.e. cxilr6716)
		--run
			Run Number
		-q|--queue
			Queue to use on SLURM
		-l|--locally
			If specified, will run locally
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
		-q|--queue)
			QUEUE="$2"
			shift
			shift
			;;
		-l|--locally)
			LOCALLY=true
			shift
			;;
		-i|--interactive)
			INTERACTIVE=true
			shift
			;;
                 *)
                        POSITIONAL+=("$1")
			shift
			;;                     
	esac
done
set -- "${POSITIONAL[@]}"

#Define cores if we don't have them
#QUEUE=${QUEUE:='psanaq'}
QUEUE=${QUEUE:='ffbh3q'}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

source /reg/g/psdm/etc/psconda.sh
ABS_PATH=/reg/g/psdm/sw/tools/smalldata_tools/examples
if [[ -v LOCALLY ]]; then
    ABS_PATH=$DIR
fi

if [[ -v INTERACTIVE ]]; then
    echo 88888
    $ABS_PATH/../examples/DataqualityPlots.py $@
else
    echo 99999
    sbatch $ABS_PATH/../examples/DataqualityPlots.py $@
fi
