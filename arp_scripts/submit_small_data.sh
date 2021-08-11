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
		-x|--norecorder
			If specified, don\'t use recorder data
		-D|--default
			If specified, translate only small data
		-f|--full
			If specified, translate everything
                -i|--image
			If specified, translate everything & save area detectors as images
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
                        EXP=$2
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


umask 002 # set permission of newly created files and dir to 664 (rwxrwxr--)

# Source the right LCLS-I/LCLS-2 stuff based on the experiment name
EXP="${EXPERIMENT:=$EXP}" # look for the environment variable for if submitted from the elog
HUTCH=${EXP:0:3}
LCLS2_HUTCHES="rix, tmo"
if echo $LCLS2_HUTCHES | grep -iw $HUTCH > /dev/null; then
    echo "This is a LCLS-II experiment"
    source /cds/sw/ds/ana/conda2/manage/bin/psconda.sh
    conda deactivate
    conda activate ps-4.2.6
    PYTHONEXE=smd2_producer.py
else
    echo "This is a LCLS-I experiment"
    source /reg/g/psdm/etc/psconda.sh -py3
    PYTHONEXE=smalldata_producer.py
fi

ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/examples/g`

#run all imports on batch node before calling mpirun on that node.
$ABS_PATH/preimport.py # do we still need that?
export PS_SRV_NODES 1 # what is that?
if [ -v NEVENTS ] && [ $NEVENTS -lt 20 ]; then
    python -u $ABS_PATH/$PYTHONEXE $@
else
    mpirun python -u $ABS_PATH/$PYTHONEXE $@
fi
