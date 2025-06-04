#!/bin/bash

usage()
{
cat << EOF
$(basename "$0"):
	Primary script to launch a smalldata_tools run analysis for LCLS-II.

	OPTIONS:
        -h|--help
            Definition of options
        -e|--experiment
            Experiment name (e.g. rixx1017523)
        -r|--run
            Run number (e.g. 520)
        --eb_cores
            Number of cores for psana EventBuilder nodes (default: 8 ; 2 for interactive)
        --bd_cores
            Number of cores for psana Big Data nodes per EB node (default: 12 ; 3 for interactive)
        --s3df
            Force S3DF data location
        --interactive
            Run interactively (default: false)
        *)
            All other arguments are passed to the cube python script.
            See the lets_cube.py help for more information.
EOF
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    case "$1" in
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
    --eb_cores)
        EB_CORES=$2
        shift
        shift
        ;;
    --bd_cores)
        BD_CORES=$2
        shift
        shift
        ;;
    --s3df)
        FORCE_S3DF=1
        shift
        ;;
    --interactive)
        INTERACTIVE=1
        shift
        ;;
    *)
        POSITIONAL+=("$1")
        shift
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# Example experiment and run numbers
    # (2mio events)
    # EXP="rixx1017523"
    # RUN="520"

    # Step scan (3mio events)
    # EXP="rixx1016923"
    # RUN="152"

    # Energy scan (3mio events, with scan PVs, no xgmd)
    # EXP="rixl1032923"
    # RUN="169"

EXP="${EXPERIMENT:=$EXP}" # default to the environment variable if submitted from the elog
RUN="${RUN_NUM:=$RUN}" # default to the environment variable if submitted from the elog
HUTCH=${EXP:0:3}
# Export EXP and RUN for when running form the CLI
# This should just re-write the existing variables when running from the elog
export EXPERIMENT=$EXP
export RUN_NUM=$RUN

# Export useful path
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
export SMD_ROOT=`echo $CWD | sed  "s|/arp_scripts||g"`

export SIT_ENV_DIR="/sdf/group/lcls/ds/ana"
echo "Sourcing environment"
source $SIT_ENV_DIR/sw/conda2/manage/bin/psconda.sh
# source /sdf/home/e/espov/dev/lcls2/setup_env.sh  # dev


# Figure out the right base path for the data (or use S3DF in force case)
if [ -v FORCE_S3DF ]; then
    DATAPATH="/sdf/data/lcls/ds"
else
    DATAPATH=`python $SMD_ROOT/arp_scripts/file_location.py -e $EXP -r $RUN`
fi
export SIT_PSDM_DATA=$DATAPATH

PYTHONEXE="lets_cube.py"
PYTHONEXE=${SMD_ROOT}/lcls2_producers/${PYTHONEXE}

if [ -v INTERACTIVE ]; then
    # Default config that should work for machine with few cores
    EB_CORES=2
    BD_CORES=3
    SRV_CORES=1
    TASKS=$(expr $BD_CORES \* $EB_CORES + $EB_CORES + $SRV_CORES + 1)  # BD + EB + SRV + SMD0

    export PS_SRV_NODES=$SRV_CORES
    export PS_EB_NODES=$EB_CORES
    
    # Run mpi locally
    MPI_CMD="mpirun -np $TASKS python -u -m mpi4py.run ${PYTHONEXE} $*"
    echo $MPI_CMD
    $MPI_CMD
    exit 0
fi


# SLURM / MPI parameters
DEFPARTITION='milano'
PARTITION=${PARTITION:=$DEFPARTITION}
ACCOUNT=${ACCOUNT:="lcls:$EXP"}

EB_CORES=${EB_CORES:=8}
BD_CORES=${BD_CORES:=12}
SRV_CORES=1
TASKS=$(expr $BD_CORES \* $EB_CORES + $EB_CORES + $SRV_CORES + 1)  # BD + EB + SRV + SMD0

export PS_SRV_NODES=$SRV_CORES
export PS_EB_NODES=$EB_CORES

# Run on slurm
LOGFILE="cube_%J.log"
SBATCH_ARGS="--use-min-nodes --exclusive -o $LOGFILE --account $ACCOUNT -p $PARTITION --ntasks $TASKS"
MPI_CMD="mpirun python -u -m mpi4py.run ${PYTHONEXE} $*"
echo "$SBATCH_ARGS" --wrap="$MPI_CMD"
sbatch $SBATCH_ARGS --wrap="$MPI_CMD"