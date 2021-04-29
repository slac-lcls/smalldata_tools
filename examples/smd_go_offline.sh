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
if [ $(echo $HOSTNAME | grep -ic -e "drp-srcf") -eq 0 ]
then
    echo "Should be run from a FFB node."
    exit
fi

EXP=${EXP:=0}
if [ $EXP -eq 0 ];
then
    echo "No experiment name given. Exit."
    exit
fi

HUTCH=${EXP:0:3}
FFB_BASE="/cds/data/drpsrcf/$HUTCH/$EXP/scratch"
PSANA_BASE="/cds/data/psdm/$HUTCH/$EXP"

# update FFB repo and get it on psana
git -C $FFB_BASE/smalldata_tools add .
git -C $FFB_BASE/smalldata_tools commit -m "Before going offline"
git -C $FFB_BASE/smalldata_tools push
# git -C $PSANA_BASE/results/smalldata_tools pull

# change jobs definition (to do)
# update_job_offline.py

# copy data over to the offline dir (probably not on us to do it)
# scp -r $FFB_BASE/hdf5/smalldata/ psana:$PSANA_BASE/hdf5/smalldata/
# scp -r $FFB_BASE/hdf5/cube/ psana:$PSANA_BASE/hdf5/smalldata/

