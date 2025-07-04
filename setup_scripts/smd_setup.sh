#!/bin/bash

usage()
{
cat << EOF
$(basename "$0"):
	Script to setup smalldata_tools on the S3DF for a given experiment.

	OPTIONS:
        -h|--help
            Definition of options
        -e
            Experiment name (i.e. cxilr6716)
        -q
            Queue. Jobs are not setup if a queue is not given
        -C|--config
            Config to use without the prod_config prefix. Default: hutch config
        --cube <config_name>
            Make cube job.
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
    -e|--experiment)
        EXP="$2"
        ARGS+=("--experiment $2")
        shift 2
        ;;
    -q|--queue)
        QUEUE="$2"
        ARGS+=("--queue $2")
        shift 2
        ;;
    -C|--config)
        CONFIG="$2"
        ARGS+=("--config $2")
        shift 2
        ;;
    --cube)
       if [[ -z "$2" ]] || [[ "$2" =~ ^- ]]; then
          echo "Error: --cube requires an argument." >&2
          exit 1
        fi
        CUBE="$2"
        ARGS+=("--cube $2")
        shift 2
        ;;
    --psplot_live)
        PSPLOT_LIVE=1
        ARGS+=("--psplot_live")
        shift
        ;;
    *)
        # all other possibilities
        ARGS+=("$1")
        shift
        ;;
    esac
done
set -- "${ARGS[@]}"

umask 002

QUEUE=${QUEUE:=0}
CUBE=${CUBE:=0}
PSPLOT_LIVE=${PSPLOT_LIVE:=0}
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" # gets the script directory.

# check that the script is run on relevant nodes
if [ $(echo $HOSTNAME | grep -ic -e "sdf") -eq 0 ]
then
    echo "This script should be run from a S3DF node. Exiting now."
    exit
fi

HUTCH=${EXP:0:3}
FFB_BASE="/sdf/data/lcls/drpsrcf/ffb/$HUTCH/$EXP"
SDF_BASE="/sdf/data/lcls/ds/$HUTCH/$EXP"

# Exit if directories dont exist
if [ ! -d "$SDF_BASE" ]; then
    echo "Experiment folder does not exist. Exiting now."
    exit
fi
if [ ! -d "$FFB_BASE" ]; then
    echo "FFB directory or mount not available for this experiment."
fi

# Clone smalldata code
echo "Cloning smalldata to S3DF experiment directory..."
if [ -d "$SDF_BASE/results/smalldata_tools" ]; then
    echo "Smalldata_tools already deployed. Skipping this step."
else
    git clone https://github.com/slac-lcls/smalldata_tools.git $SDF_BASE/results/smalldata_tools
fi
echo "... Done."

# change smalldata_tool permissions for the ARP (temporary?)
# ARP kubernetes pods do not have the ACL
chmod -R o+r $SDF_BASE/results/smalldata_tools
chmod o+x $SDF_BASE/results/smalldata_tools
chmod -R o+x $SDF_BASE/results/smalldata_tools/arp_scripts
chmod -R o+x $SDF_BASE/results/smalldata_tools/lcls1_producers
chmod -R o+x $SDF_BASE/results/smalldata_tools/lcls2_producers

# Create h5 and plot directories
mkdir -p $SDF_BASE/hdf5/smalldata
mkdir -p $SDF_BASE/hdf5/smalldata/cube
mkdir -p $SDF_BASE/stats/summary/Cube

# make arp jobs
echo "Now making ARP jobs"
LCLS1_HUTCHES=("xcs" "cxi" "mec")
LCLS2_HUTCHES=("tmo" "txi" "rix" "ued" "mfx" "xpp")

if [ $QUEUE != "0" ]; then
    if [[ ${LCLS1_HUTCHES[@]} =~ $HUTCH ]]; then
        source /sdf/group/lcls/ds/ana/sw/conda1/manage/bin/psconda.sh
        python $MYDIR/make_arp_jobs_lcls1.py $*
    elif [[ ${LCLS2_HUTCHES[@]} =~ $HUTCH ]]; then
        source /sdf/group/lcls/ds/ana/sw/conda2/manage/bin/psconda.sh
        python $MYDIR/make_arp_jobs_lcls2.py $*
    fi
fi
