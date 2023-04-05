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
			Number of events to analyze
		-q|--queue
			Queue to use on SLURM
		-c|--cores
			Number of cores to be utilized
		-f|--full
			If specified, translate everything
		-D|--default
			If specified, translate only smalldata
                -i|--image
			If specified, translate everything & save area detectors as images
		--norecorder
			If specified, don't use recorder data
                --nparallel
                        Number of processes per node
                --postTrigger
                        Post that primary processing done to elog to seconndary jobs can start
                --interactive
                        Run the process live w/o batch system
                --logdir
                        save log-files in specified directory
EOF

}

#Look for args we expect, for now ignore other args
#Since we can have a mix of flags/args and args do in loop

#Generate smalldata_producer_arp.py args, for now if
#Experiment and Run not specified assume it's being handed
#by the ARP args in the os.environment

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
    --nparallel)
        TASKS_PER_NODE="$2"
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
        POSITIONAL+=("--experiment $2")
        EXP=$2
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
    *)
        POSITIONAL+=("$1")
        shift
        ;;                     
    esac
done
set -- "${POSITIONAL[@]}"


umask 002 # set permission of newly created files and dir to 664 (rwxrwxr--)

# Source the right LCLS-I/LCLS-2 stuff based on the experiment name
EXP="${EXPERIMENT:=$EXP}" # default to the environment variable if submitted from the elog
HUTCH=${EXP:0:3}
LCLS2_HUTCHES="rix, tmo, ued"
SIT_ENV_DIR="/cds/group/psdm"
S3DF="sdf"

export MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/producers/g`


DEFQUEUE='psanaq'
if [[ $HOSTNAME == *drp* ]]; then
    DEFQUEUE='anaq'
elif [[ $HOSTNAME == *sdf* ]]; then
    DEFQUEUE='milano'
fi
#Define cores if we don't have them
#Set to 1 by default
CORES=${CORES:=1}
QUEUE=${QUEUE:=$DEFQUEUE}
ACCOUNT=${ACCOUNT:='lcls'}
RESERVATION=${RESERVATION:=''}
#QUEUE=${QUEUE:='anaq'}
# select tasks per node to match the number of cores:
if [[ $QUEUE == *psanaq* ]]; then
    TASKS_PER_NODE=${TASKS_PER_NODE:=12}
elif [[ $QUEUE == *psfeh* ]]; then
    TASKS_PER_NODE=${TASKS_PER_NODE:=16}
elif [[ $QUEUE == *ffb* ]]; then
    TASKS_PER_NODE=${TASKS_PER_NODE:=60}
elif [[ $QUEUE == *milano* ]]; then
    TASKS_PER_NODE=${TASKS_PER_NODE:=128}
else
    TASKS_PER_NODE=${TASKS_PER_NODE:=12}
fi

if [ $TASKS_PER_NODE -gt $CORES ]; then
    TASKS_PER_NODE=$CORES 
fi

if [[ $HOSTNAME == *sdf* ]]; then
    SIT_ENV_DIR="/sdf/group/lcls/ds/ana"
fi
if echo $LCLS2_HUTCHES | grep -iw $HUTCH > /dev/null; then
    echo "This is a LCLS-II experiment"
    source $SIT_ENV_DIR/sw/conda2/manage/bin/psconda.sh
    PYTHONEXE=smd2_producer.py
    export PS_SRV_NODES=1 # 1 is plenty enough for the 120 Hz operation
else
    echo "This is a LCLS-I experiment"
    #echo "Setting up the enviroment: "$SIT_ENV_DIR/sw/ds/ana/conda1/manage/bin/psconda.sh
    if [[ $HOSTNAME == *sdf* ]]; then
        source $SIT_ENV_DIR/sw/conda1/manage/bin/psconda.sh
    else
        source /cds/sw/ds/ana/conda1/manage/bin/psconda.sh
    fi
    PYTHONEXE=smd_producer.py
fi

echo ---- print environment ----
env | sort
echo --- printed environment ---

if [ -v INTERACTIVE ]; then
    #run all imports on batch node before calling mpirun on that node.
    if [ -v NEVENTS ] && [ $NEVENTS -lt 20 ]; then
        python -u $ABS_PATH/$PYTHONEXE $@
    else
        mpirun -np $CORES python -u $ABS_PATH/$PYTHONEXE $@
    fi

    exit 0
fi

LOGFILE='smd_'${EXP}'_Run'${RUN}'_%J.log'
if [ -v LOGDIR ]; then
    if [ ! -d "$LOGDIR" ]; then
        mkdir -p "$LOGDIR"
    fi
    LOGFILE=$LOGDIR'/'$LOGFILE
fi

if [[ $QUEUE == *milano* ]]; then
    if [[ $ACCOUNT == 'lcls' ]]; then
	sbatch -p $QUEUE --ntasks-per-node $TASKS_PER_NODE --ntasks $CORES --exclusive --account $ACCOUNT --qos preemptable -o $LOGFILE --wrap="mpirun -np $CORES python -u ${ABS_PATH}/${PYTHONEXE} $*"
    else
        echo ---- $ABS_PATH/$PYTHONEXE $@
	sbatch -p $QUEUE --ntasks-per-node $TASKS_PER_NODE --ntasks $CORES --exclusive --account $ACCOUNT -o $LOGFILE --wrap="mpirun -np $CORES python -u ${ABS_PATH}/${PYTHONEXE} $*"

    fi
else
    sbatch -p $QUEUE --ntasks-per-node $TASKS_PER_NODE --ntasks $CORES --exclusive -o $LOGFILE --wrap "mpirun -np $CORES python -u ${ABS_PATH}/${PYTHONEXE} $*"
fi
