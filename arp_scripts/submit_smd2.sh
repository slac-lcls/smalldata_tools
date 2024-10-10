#!/bin/bash

usage()
{
cat << EOF
$(basename "$0"):
	Primary script to launch a smalldata_tools run analysis for LCLS-II.

	OPTIONS:
        -h|--help
            Definition of options
EOF
}

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
        EXP=$2
        POSITIONAL+=("--experiment $2")
        shift
        shift
        ;;
    -r|--run)
        RUN=$2
        POSITIONAL+=("--run $2")
        shift
        shift
        ;;
    -d|--directory)
        POSITIONAL+=("--directory $2")
        shift
        shift
        ;;
    -n|--nevents)
        NEVENTS=$2
        POSITIONAL+=("--nevents $2")
        shift
        shift
        ;;
    --interactive)
        INTERACTIVE=1
        shift
        ;;
    --nodes)
        NODES=$2
        shift
        shift
        ;;
    --eb_cores)
        EB_CORES=$2
        shift
        shift
        ;;
    --srv_cores)
        SRV_CORES=$2
        shift
        shift
        ;;
    --account)
        ACCOUNT="$2"
        shift
        shift
        ;;
    --reservation)
        RESERVATION="$2"
        shift
        shift
        ;;
    --logdir)
        LOGDIR="$2"
        shift
        shift
        ;;
    -p|--partition)
        PARTITION="$2"
        shift
        shift
        ;;
    --s3df)
        FORCE_S3DF=1
        shift
        ;;
    --mpi_optim)
        MPI_OPTIM=1
        POSITIONAL+=("--mpi_optim")
        shift
        ;;
    --psplot_live)
        PSPLOT_LIVE=1
        POSITIONAL+=("--psplot_live_mode")
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

EXP="${EXPERIMENT:=$EXP}" # default to the environment variable if submitted from the elog
RUN="${RUN_NUM:=$RUN}" # default to the environment variable if submitted from the elog
HUTCH=${EXP:0:3}
ARP_LOCATION="${ARP_LOCATION:=LOCAL}"

# Export EXP and RUN for when running form the CLI
# This should just re-write the existing variables when running from the elog
export EXPERIMENT=$EXP
export RUN_NUM=$RUN

# Export useful path
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
export SMD_ROOT=`echo $MYDIR | sed  "s|/arp_scripts||g"`

export SIT_ENV_DIR="/sdf/group/lcls/ds/ana"

# Source env. Needed to get python
echo "Sourcing LCLS-II environment"
source $SIT_ENV_DIR/sw/conda2/manage/bin/psconda.sh
# source /sdf/home/e/espov/dev/lcls2/setup_env.sh

# Figure out the right base path for the data (or use S3DF in force case)
if [ -v FORCE_S3DF ]; then
    DATAPATH="/sdf/data/lcls/ds"
else
    DATAPATH=`python $SMD_ROOT/arp_scripts/file_location.py -e $EXP -r $RUN`
fi
export SIT_PSDM_DATA=$DATAPATH

if [ -v INTERACTIVE ]; then
    export PS_SRV_NODES=1
    export PS_EB_NODES=1
    $SMD_ROOT/arp_scripts/run_smd2.sh $*
    exit 0
fi

# SLURM / MPI parameters
DEFPARTITION='milano'
PARTITION=${PARTITION:=$DEFPARTITION}
ACCOUNT=${ACCOUNT:="lcls:$EXP"}


# EB, BD and SRV core allocation
if [[ "$PARTITION" == "milano" ]]; then
    # For milano:
    CORES_PER_NODE=120
    DEFAULT_NODES=2
    NODES=${NODES:=$DEFAULT_NODES}

    if [ -v MPI_OPTIM ]; then
        MPI_SLOTS=$(($CORES_PER_NODE*($NODES-1)))  # Number of cores on all nodes but the one for SMD0
        DEFAULT_SRV_CORES=$((16*($NODES-1)))  # 16 writers per node minus the SMD0 node
    else
        MPI_SLOTS=$(($CORES_PER_NODE*$NODES-1))  # Total number of cores minus SMD0
        DEFAULT_SRV_CORES=$((16*$NODES))  # 16 writers per node seems to be a good number for milano
    fi

    if [ -v PSPLOT_LIVE ]; then
        # Use a single SRV node for psplot_live
        DEFAULT_SRV_CORES=1
    fi

    SRV_CORES=${SRV_CORES:=$DEFAULT_SRV_CORES}

    DEFAULT_EB_CORES=$((($MPI_SLOTS-$SRV_CORES)/16))  # 1/16 is the recommended ratio for generic jobs.
    EB_CORES=${EB_NODES:=$DEFAULT_EB_CORES}
fi

export PS_SRV_NODES=$SRV_CORES
export PS_EB_NODES=$EB_CORES

SBATCH_ARGS="--nodes $NODES --account $ACCOUNT -p $PARTITION"

echo "sbatch arguments: $SBATCH_ARGS"

sbatch $SBATCH_ARGS $SMD_ROOT/arp_scripts/run_smd2.sh $*
