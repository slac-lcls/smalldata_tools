#!/bin/bash

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
        -n|--nevents
            Number of events per bin to analyze
        -q|--queue
            Queue to use on SLURM
        -c|--cores
            Number of cores to be utilized
        --maxnodes
            Max number of nodes to use
        --interactive
            Run the process live w/o batch system
        --logdir
            save log-files in specified directory
        --s3df
            Forces to load xtc files from the S3DF location
        --reservation
            Slurm reservation name
EOF

}
#Look for args we expect, for now ignore other args
#Since we can have a mix of flags/args and args do in loop

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
        usage
        exit
        ;;
    -q|--queue)
        QUEUE="$2"
        shift
        shift
        ;;
    -c|--cores)
        CORES="$2"
        shift
        shift
        ;;
    --maxnodes)
        MAX_NODES="$2"
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
    --interactive)
        INTERACTIVE=1
        shift
        ;;
    --s3df)
        FORCE_S3DF=1
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
RUN="${RUN_NUM:=$RUN}" # same as EXP
HUTCH=${EXP:0:3}
SIT_ENV_DIR="/cds/group/psdm"
ARP_LOCATION="${ARP_LOCATION:=LOCAL}"

# Export EXP and RUN for when running form the CLI
# This should just re-write the existing variables when running from the elog
export EXPERIMENT=$EXP
export RUN_NUM=$RUN

export MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/lcls1_producers/g`

# Default queue and S3DF flag
DEFQUEUE='psanaq'
if [[ $HOSTNAME == *drp* ]]; then
    DEFQUEUE='anaq'
fi

if [ -d "/sdf/data/lcls/" ]; then
    ON_S3DF=true
    SIT_ENV_DIR="/sdf/group/lcls/ds/ana"
    DEFQUEUE='milano'
else
    ON_S3DF=false
fi

# Define # of cores
QUEUE=${QUEUE:=$DEFQUEUE}
ACCOUNT=${ACCOUNT:="lcls:$EXP"}

# Select number of cores and max nodes
if [[ $QUEUE == *milano*  ]]; then
    CORES=${CORES:=120} # a full node by default
    MAX_NODES=${MAX_NODES:=4}
else
    CORES=${CORES:=1} # default to 1 outside S3DF
    MAX_NODES=${MAX_NODES:=4}
fi

# Source the right LCLS-I/LCLS-2 stuff based on the experiment name
echo "This is a LCLS-I experiment"
if $ON_S3DF; then
    source $SIT_ENV_DIR/sw/conda1/manage/bin/psconda.sh
else
    source /cds/sw/ds/ana/conda1/manage/bin/psconda.sh
fi
PYTHONEXE=letsCube.py

# figure out the right base path for the data (or use S3DF in force case)
if [ -v FORCE_S3DF ]; then
    DATAPATH="/sdf/data/lcls/ds"
else
    DATAPATH=`python $MYDIR/file_location.py -e $EXP -r $RUN`
fi
export SIT_PSDM_DATA=$DATAPATH

echo "SIT_PSDM_DATA: $SIT_PSDM_DATA"


if [ -v INTERACTIVE ]; then
    # run in local terminal
    CORES=8 # use only 8 cores in interactive mode
    mpirun -np $CORES python -u -m mpi4py.run $ABS_PATH/$PYTHONEXE $@
    exit 0
fi

LOGFILE='cube_'${EXPERIMENT}'_Run'${RUN_NUM}'_%J.log'
if [ -v LOGDIR ]; then
    if [ ! -d "$LOGDIR" ]; then
        mkdir -p "$LOGDIR"
    fi
    LOGFILE=$LOGDIR'/'$LOGFILE
fi

SBATCH_ARGS="-p $QUEUE --nodes 0-$MAX_NODES --ntasks $CORES -o $LOGFILE"
MPI_CMD="mpirun -np $CORES python -u -m mpi4py.run ${ABS_PATH}/${PYTHONEXE} $*"


if [[ $QUEUE == *milano* ]]; then
    if [ -v RESERVATION ]; then
        SBATCH_ARGS="$SBATCH_ARGS --reservation $RESERVATION"
    fi
    if [[ $ACCOUNT == 'lcls' ]]; then
        echo $SBATCH_ARGS --qos preemtable --account $ACCOUNT --wrap="$MPI_CMD"
        sbatch $SBATCH_ARGS --qos preemptable --account $ACCOUNT --wrap="$MPI_CMD"
    else
        echo ---- $ABS_PATH/$PYTHONEXE $@
        echo $SBATCH_ARGS --account $ACCOUNT --wrap="$MPI_CMD"
        sbatch $SBATCH_ARGS --account $ACCOUNT --wrap="$MPI_CMD"

    fi
else # for outside s3df
    sbatch $SBATCH_ARGS --wrap "$MPI_CMD"
fi

