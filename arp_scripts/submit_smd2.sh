#!/bin/bash


usage()
{
cat << EOF
$(basename "$0"):
    Primary script to launch a smalldata_tools run analysis for LCLS-II.

    OPTIONS:
        -h, --help
            Show this help message and exit.

        -e, --experiment EXP
            Specify the experiment name (e.g., cxilr6716).

        -r, --run RUN
            Specify the run number.

        -d, --directory DIR
            Specify the output directory.

        -n, --nevents NEVENTS
            Number of events to process.

        --interactive
            Run in interactive mode (no SLURM submission).

        --nodes NODES
            Number of nodes to use (default: 2 for milano).

        --ntasks NTASKS
            Total number of MPI tasks to use.

        --memory MEMORY
            Memory allocation per task. Only used if --ntasks is specified.

        --eb_cores EB_CORES
            Number of event builder cores.

        --srv_cores SRV_CORES
            Number of server cores.

        --account ACCOUNT
            SLURM account to use (default: lcls:\$EXP).

        --reservation RESERVATION
            SLURM reservation to use.

        --logdir LOGDIR
            Directory for log files.

        -p, --partition PARTITION
            SLURM partition to use (default: milano).

        --s3df
            Force use of S3DF data path.

        --mpi_optim
            Enable MPI optimization for very high rates.

        --psplot_live
            Enable psplot live mode (single SRV node).

    Example usage:
        $(basename "$0") -e cxilr6716 -r 45 --nodes 4

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
    --ntasks)
        NTASKS=$2
        shift
        shift
        ;;
    --memory)
        MEMORY=$2
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


# If ntasks is specified, nodes and srv_cores must be specified too. This
#  is the fully custom mode.
if [[ -v NTASKS && (! -v NODES || ! -v SRV_CORES) ]]; then
  echo "Error: --ntasks requires both --nodes and --srv_cores to be also specified."
  exit 1
fi


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
# source /sdf/group/lcls/ds/ana/sw/conda2/manage/bin/pscondatest.sh  # test env for new interface

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

SBATCH_ARGS="--account $ACCOUNT -p $PARTITION"

# EB, BD and SRV core allocation
# Full custom
if [ -v NTASKS ] && [ -v NODES ]; then
    echo "Using user-defined NTASKS: $NTASKS, NODES: $NODES, SRV_CORES: $SRV_CORES"
    MPI_SLOTS=$NTASKS

    SBATCH_ARGS+=" --nodes=$NODES --ntasks=$NTASKS"
    if [ -v MEMORY ]; then
        SBATCH_ARGS+=" --mem-per-cpu=$MEMORY"
    fi

# For milano:
elif [[ "$PARTITION" == "milano" ]]; then
    CORES_PER_NODE=120
    DEFAULT_NODES=2
    NODES=${NODES:=$DEFAULT_NODES}

    SBATCH_ARGS+=" --nodes=$NODES"

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
fi

SRV_CORES=${SRV_CORES:=$DEFAULT_SRV_CORES}

DEFAULT_EB_CORES=$((($MPI_SLOTS-$SRV_CORES)/16))  # 1/16 is the recommended ratio for generic jobs.
EB_CORES=${EB_NODES:=$DEFAULT_EB_CORES}

export PS_SRV_NODES=$SRV_CORES
export PS_EB_NODES=$EB_CORES

# If we run into the chunk size error:
# export PS_SMD_CHUNKSIZE=1073741824

echo "sbatch arguments: $SBATCH_ARGS"

sbatch $SBATCH_ARGS $SMD_ROOT/arp_scripts/run_smd2.sh $*
