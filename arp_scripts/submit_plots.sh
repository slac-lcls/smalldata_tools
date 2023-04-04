#!/bin/bash

# Took this from Silke's script, nice implementation
usage()
{
cat << EOF
$(basename "$0"): 
	Script to launch a smalldata_tools run analysis
	
	OPTIONS:
		-h|--help
			Definition of options
		--experiment
			Experiment name (i.e. cxilr6716)
		--run
			Run Number
		-q|--queue
			Queue to use on SLURM
		-i|--interactive
			Run interactively
		-p|--pedestals
			plot the pedestals instead
		-b|--bld
			output selected BLD data in a text file (MEC style)
EOF

}
#Look for args we expect, for now ignore other args
#Since we can have a mix of flags/args and args do in loop

#Generate smalldata_producer_arp.py args, for now if
#Experiment and Run not specified assume it's being handed
#by the ARP args in the os.environment

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
		-p|--pedestals)
			PEDESTAL=true
			shift
			;;
		-b|--bld)
			BLD=true
			shift
			;;
        *)
            POSITIONAL+=("$1")
			shift
			;;                     
	esac
done
set -- "${POSITIONAL[@]}"

#Define cores if we don't have them
QUEUE=${QUEUE:='psanaq'}
#QUEUE=${QUEUE:='ffbh3q'}

export MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

EXP="${EXPERIMENT:=$EXP}" # default to the environment variable if submitted from the elog
HUTCH=${EXP:0:3}
LCLS2_HUTCHES="rix, tmo, ued"
if echo $LCLS2_HUTCHES | grep -iw $HUTCH > /dev/null; then
    echo "This is a LCLS-II experiment"
    source /cds/sw/ds/ana/conda2/manage/bin/psconda.sh
    ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/summaries/g`
    PLOT_PY=BeamlineSummaryPlots_$HUTCH
    if [[ -v PEDESTAL ]]; then
        PLOT_PY=PedestalPlot
    fi
else
    source /reg/g/psdm/etc/psconda.sh -py3
    #conda activate ana-4.0.16-py3
    ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/summaries/g`
    PLOT_PY=BeamlineSummaryPlots_$HUTCH
    if [[ -v PEDESTAL ]]; then
        PLOT_PY=PedestalPlot
    elif [[ -v BLD ]]; then
        ABS_PATH=`echo $MYDIR | sed  s/arp_scripts/producers/g`
        PLOT_PY=BldEpics
    fi
fi

echo calling $ABS_PATH/$PLOT_PY.py $@
if [[ -v INTERACTIVE ]]; then
    $ABS_PATH/$PLOT_PY.py $@
else
    sbatch -p $QUEUE $ABS_PATH/$PLOT_PY.py $@
fi
