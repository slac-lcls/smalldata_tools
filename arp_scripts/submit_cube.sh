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
		-q|--queue
			Queue to use on SLURM
		-c|--cores
			Number of cores to be utilized
		-d|--directory
			Full path to directory for output file
		-n|--nevents
			Number of events to analyze
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

#
#use the old analysis release to set up with old mpi-parallel hdf5 writing.
#
source /reg/g/psdm/etc/psconda.sh.old
ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/examples/g`

if [ -v NEVENTS ] && [ $NEVENTS -lt 10 ]; then
    $ABS_PATH/MakeCube.py $@
else
    mpirun $ABS_PATH/MakeCube.py $@
fi
