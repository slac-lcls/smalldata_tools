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
		---directory
			Full path to directory for output file
		--nevents
			Number of events to analyze
		--norecorder
			If specified, don't use recorder data
		-q|--queue
			Queue to use on SLURM
		-c|--cores
			Number of cores to be utilized
		-s|--single
			Run on a single core
		-l|--locally
			If specified, will run locally
		-F|--full
			eIf specified, translate everything
		-t|--test
			Run the slurm job as test only to get job info
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
                 *)
                        POSITIONAL+=("$1")
			shift
			;;                     
	esac
done
set -- "${POSITIONAL[@]}"

#Define cores if we don't have them
#Set to 1 if single is set
CORES=${CORES:=1}
QUEUE=${QUEUE:='psanaq'}

# deal with request cores & cores/per nodde
#if [ $CORES -gt 1 ]; then
#   if [ $QUEUE =~ 'psneh' ]; then

#SBATCH --ntasks-per-node=$CORES
#SBATCH --nodes=1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
sbatch -p $QUEUE --ntasks-per-node $CORES $DIR/submit_small_data_inner.sh $@
