#!/bin/bash

usage()
{
cat << EOF
$(basename "$0"): 
	Script to setup smalldata_tools for an experiment.
	
	OPTIONS:
		-h|--help
			Definition of options
		-e
			Experiment name (i.e. cxilr6716)
        -q
            Queue (only for ffb)
        --nopsana
            Do not setup smalldata on psana.
        --noffb
        	Do not setup smalldata on the FFB.
            If exists on psana, will clone the repo from psana.
EOF
}

ARGS=()
while [[ $# -gt 0 ]]
do
    KEY="$1"
    case $KEY in
    -h|--help)
		usage
		exit
		;;
    -e)
        EXP="$2"
        shift 2
        ;;
    -q)
        QUEUE="$2"
        shift 2
        ;;
    --noffb)
        FFB=0
        shift 1
        ;;
    --nopsana)
        PSANA=0
        shift 1
        ;;
    *) # all other possibilities
        ARGS+=("$1")
        echo $@
        shift
        ;;
    esac
done

FFB=${FFB:=1}
PSANA=${PSANA:=1}
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" # gets the script directory.

# check that the script is run on relevant nodes
if [ $FFB -eq 1 ]; then
    if [ $(echo $HOSTNAME | grep -ic -e "drp-srcf") -eq 0 ]
    then
        echo "Should be run from a FFB node."
        exit
    fi
elif [ $PSANA -eq 1 ]; then
    if [ $(echo $HOSTNAME | grep -ic -e "drp-srcf" -e "psana") -eq 0 ]
    then
        echo "Should be run from a FFB or PSANA node."
        exit
    fi
fi

HUTCH=${EXP:0:3}
FFB_BASE="/cds/data/drpsrcf/$HUTCH/$EXP/scratch"
PSANA_BASE="/cds/data/psdm/$HUTCH/$EXP"

# Exit if directories dont exist
if [ $FFB -eq 1 ]; then
    if [ ! -d "$FFB_BASE" ]; then
      exit
    fi
fi
if [ $PSANA -eq 1 ]; then
    if [ ! -d "$PSANA_BASE" ]; then
      exit
    fi
fi

# setup smalldata code
# On PSANA
if [ $PSANA -eq 1 ]; then
    echo "Cloning smalldata to PSANA experiment directory..."
    if [ -d "$PSANA_BASE/results/smalldata_tools" ]; then
        echo "Smalldata_tools already on PSANA"
    else
        git clone https://github.com/slac-lcls/smalldata_tools.git $PSANA_BASE/results/smalldata_tools
    fi
    echo "... Done."
fi
echo "Sleep 2"
sleep 2
# On FFB
if [ $FFB -eq 1 ]; then
    echo "Cloning smalldata to FFB experiment directory..."
    if [ -d "$FFB_BASE/smalldata_tools" ]; then
        echo "Smalldata_tools already on FFB"
    else
        if [ -d "$PSANA_BASE/results/smalldata_tools" ]; then
            echo "Cloning from the PSANA directory."
            git clone $PSANA_BASE/results/smalldata_tools $FFB_BASE/smalldata_tools
        else
            echo "Cloning from the remote."
            git clone https://github.com/slac-lcls/smalldata_tools.git $FFB_BASE/smalldata_tools
        fi
    fi
    echo "... Done."
fi


# Create h5 directories
if [ $FFB -eq 1 ]; then
    mkdir -p $FFB_BASE/hdf5/smalldata
    mkdir -p $FFB_BASE/hdf5/cube
    mkdir -p $FFB_BASE/hdf5/debug
fi
# if [ $PSANA -eq 1 ]; then
#     mkdir -p $PSANA_BASE/hdf5/smalldata/debug
#     mkdir -p $PSANA_BASE/hdf5/cube
# fi

# make arp jobs
python $MYDIR/make_arp_jobs.py --experiment $EXP --queue $QUEUE --psana $PSANA --ffb $FFB
