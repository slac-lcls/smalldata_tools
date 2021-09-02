#!/bin/bash

usage()
{
cat << EOF
$(basename "$0"): 
	Script to move data and data analysis code to the offline (psana) directories
	
	OPTIONS:
		-h|--help
			Definition of options
		-e
			Experiment name (i.e. cxilr6716)
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
    *) # all other possibilities
        ARGS+=("$1")
        echo $@
        shift
        ;;
    esac
done

# check that the script is run on relevant nodes
#if [ $(echo $HOSTNAME | grep -ic -e "drp-srcf") -eq 0 ]
#then
#    echo "Should be run from a FFB node."
#    exit
#fi

EXP=${EXP:=0}
if [ $EXP -eq 0 ];
then
    echo "No experiment name given. Exit."
    exit
fi

HUTCH=${EXP:0:3}
FFB_BASE="/reg/data/drpsrcf/$HUTCH/$EXP/scratch"
PSANA_BASE="/cds/data/psdm/$HUTCH/$EXP"
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" # gets the script directory.

# update FFB repo and get it on psana
git -C $PSANA_BASE/results/smalldata_tools config receive.denyCurrentBranch=updateInstead # necessary to be able to push ffb repo
git -C $FFB_BASE/smalldata_tools add .
git -C $FFB_BASE/smalldata_tools commit -m "Going offline"
git -C $FFB_BASE/smalldata_tools push -f

# change jobs definition
python $MYDIR/update_arp_jobs.py --experiment $EXP


# copy h5 files to offline
#cp -uv $FFB_BASE/hdf5/smalldata/*.h5 $PSANA_BASE/hdf5/smalldata
