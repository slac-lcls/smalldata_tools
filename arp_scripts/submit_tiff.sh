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
		-q|--queue
			Queue to use on SLURM
		-i|--interactive
			Run interactively
EOF

}

POSITIONAL=()
EXP=$EXPERIMENT
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
		-r|--run)
		    POSITIONAL+=("--run $2")
		    RUN=$2
		    shift
		    shift
		    ;;
		-e|--experiment)
		    POSITIONAL+=("--experiment $2")
		    EXP=$2
		    shift
		    shift
		    ;;
		-i|--interactive)
		    INTERACTIVE=true
		    shift
		    ;;
        *)
            POSITIONAL+=("$1")
			shift
			;;                     
	esac
done
set -- "${POSITIONAL[@]}"

SIT_ENV_DIR="/sdf/group/lcls/ds/ana"

#Define cores if we don't have them
QUEUE=${QUEUE:='milano'}
ACCOUNT=${ACCOUNT:="lcls:$EXP"}

export MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

RUN="${RUN_NUM:=$RUN}" # default to the environment variable if submitted from the elog
EXP="${EXPERIMENT:=$EXP}" # default to the environment variable if submitted from the elog
HUTCH=${EXP:0:3}

source $SIT_ENV_DIR/sw/conda1/manage/bin/psconda.sh
#conda activate ana-4.0.16-py3
ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/summaries/g`
EXE_PY=MakeTiff

echo calling $ABS_PATH/$EXE_PY.py $@
if [[ -v INTERACTIVE ]]; then
    $ABS_PATH/$EXE_PY.py $@
else
    sbatch -p $QUEUE --account $ACCOUNT $ABS_PATH/$EXE_PY.py $@
fi
