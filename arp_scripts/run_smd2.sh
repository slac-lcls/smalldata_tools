#!/bin/bash

#SBATCH --job-name=smd
#SBATCH --exclusive
# Note: when run trough sbatch, this script runs from the slurm scheduler folder

usage()
{
cat << EOF
$(basename "$0"):
	Primary script to launch a smalldata_tools run analysis for LCLS-II.

	OPTIONS:
        -h|--help
            Definition of options

        --mpi-optim
            Use special MPI setup to optimize SMD0 speed. Require to use
            nodes exclusively, and may thus lead to longer delay until resources
            are allocated

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
    --mpi_optim)
        MPI_OPTIM=1
        shift
        ;;
    *)
        POSITIONAL+=("$1")
        shift
        ;;
    esac
done
set -- "${POSITIONAL[@]}"

echo "SIT_ENV_DIR: $SIT_ENV_DIR"
echo "SIT_PSDM_DATA: $SIT_PSDM_DATA"

# echo "Sourcing LCLS-II environment"
# source $SIT_ENV_DIR/sw/conda2/manage/bin/psconda.sh
PYTHONEXE=smd_producer.py

# MPI Setup
echo "PS_EB_NODES: $PS_EB_NODES"
echo "PS_SRV_NODES: $PS_SRV_NODES"

if [ -v MPI_OPTIM ]; then
    # Optimize psana2 parallelization
    # This needs to run inside a batch job
    if [ -n $SLURM_JOB_ID ]; then
        echo "Special MPI setup"
        source $SMD_ROOT/arp_scripts/setup_hosts_openmpi.sh
    fi
fi

echo "Producer command: $SMD_ROOT/lcls2_producers/$PYTHONEXE $@"
mpirun python -u -m mpi4py.run $SMD_ROOT/lcls2_producers/$PYTHONEXE $@
