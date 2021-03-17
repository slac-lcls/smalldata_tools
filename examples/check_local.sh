#!/bin/bash

usage()
{
cat << EOF
$(basename "$0"): 
    Script to locally check the smd translation on a reduced number of events. Needs to have a /tests/ directory in the experiment hdf5 folder.
	
	OPTIONS:
		-h
			Definition of options
		-e
			Experiment name (i.e. cxilr6716)
		-r
			Run Number
        -n
            Number of events (default 50)
EOF
}

# parse arguments
while getopts hr:e:n: option # if ':' it must take an argument, 
do
    case "${option}" in 
        h)
			usage
			exit
			;;
        r) RUN=${OPTARG};; 
        e) EXP=${OPTARG};;
        n) $NEVENT=${OPTARG};;
    esac 
done

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
RUN=${RUN:=0}
if [[ $RUN == 0 ]]; then
    echo "Please give a run number (-r)"
    exit
fi
EXP=${EXP:='xpplv7918'} # change default exp to avoid having to type it all the time
NEVENTS=${NEVENT:=10}

if [[ $MYDIR == *"drpsrcf"* ]]; then
    echo "On drp"
    H5PATH="$( cd $MYDIR && cd ../../hdf5/debug >/dev/null && pwd )"
#     if [ ! -d "$H5PATH" ]; then
#         mkdir $H5PATH
#     fi
elif [[ $MYDIR == *"psdm"* ]]; then
    echo "On psana"
    H5PATH="$( cd $MYDIR && cd ../../../hdf5/debug >/dev/null && pwd )"
else
    echo "On some unknown place."
    H5PATH=$MYDIR  
fi

# echo $MYDIR
# echo $H5PATH
# echo $RUN
# echo $EXP

# /cds/home/e/espov/lcls_software_tools/smalldata_tools/arp_scripts/submit_smd.sh --interactive --nevents $NEVENTS --directory $H5PATH --experiment $EXP --run $RUN
/cds/home/e/espov/lcls_software_tools/smalldata_tools/arp_scripts/submit_smd.sh --interactive --nevents $NEVENTS --directory /cds/data/psdm/xpp/xpplv7918/hdf5/debug --experiment $EXP --run $RUN