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
        -c|--copy
            If given, will make procserv to copy h5 files from the ffb to anafs
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
    -c|--copy)
        COPY=1
        shift
        ;;
    *) # all other possibilities
        ARGS+=("$1")
        echo $@
        shift
        ;;
    esac
done

# check that the script is run on relevant nodes
if [ $(echo $HOSTNAME | grep -ic -e "psana") -eq 0 ]
then
   echo "Should be run from a psana node."
   exit
fi

EXP=${EXP:=0}
if [ $EXP -eq 0 ];
then
    echo "No experiment name given. Exit."
    exit
fi
COPY=${COPY:=0}

export HUTCH=${EXP:0:3}
export FFB_BASE="/reg/data/drpsrcf/$HUTCH/$EXP/scratch"
export PSANA_BASE="/cds/data/psdm/$HUTCH/$EXP"
export MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" # gets the script directory.

# update FFB repo and get it on psana
git -C $PSANA_BASE/results/smalldata_tools config receive.denyCurrentBranch=updateInstead # necessary to be able to push ffb repo
git -C $FFB_BASE/smalldata_tools add .
git -C $FFB_BASE/smalldata_tools commit -m "Going offline"
git -C $FFB_BASE/smalldata_tools push -f

# change jobs definition
source /reg/g/psdm/etc/psconda.sh -py3
python $MYDIR/update_arp_jobs.py --experiment $EXP


# copy h5 files to offline in a procserver
if [[ $COPY -ne 0 ]]; then
    #cp -uv $FFB_BASE/hdf5/smalldata/*.h5 $PSANA_BASE/hdf5/smalldata
    export NAME="smd_ffb_anafs_transfer"
    export PROCSERV="/cds/group/pcds/pkg_mgr/release/procServ/2.8.0-1.3.0/rhel7-x86_64/bin/procServ --oneshot --ignore ^D^C"
    #$PROCSERV --logfile $PSANA_BASE/scratch/$NAME.log --name $NAME 43400 rsync -avu $FFB_BASE/hdf5/smalldata/*.h5 $PSANA_BASE/hdf5/smalldata
    $PROCSERV --logfile $PSANA_BASE/scratch/$NAME.log --name $NAME 43400 $MYDIR/sync_h5.cmd
fi
