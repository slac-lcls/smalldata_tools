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

        -n|--cores
            Number of MPI ranks to launch with mpirun

        --gather_interval
            Number of events per smalldata gather (default: 100)

        --mpi_optim
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
    --mpi_optim|--mpi-optim)
        MPI_OPTIM=1
        shift
        ;;
    -n|--cores)
        MPI_CORES="$2"
        shift
        shift
        ;;
    --gather_interval)
        GATHER_INTERVAL="$2"
        POSITIONAL+=("--gather_interval" "$2")
        shift
        shift
        ;;
    *)
        POSITIONAL+=("$1")
        shift
        ;;
    esac
done
if [ -z "${GATHER_INTERVAL}" ]; then
    GATHER_INTERVAL=100
    POSITIONAL+=("--gather_interval" "$GATHER_INTERVAL")
fi

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
    if [ -n "$SLURM_JOB_ID" ]; then
        echo "Special MPI setup"
        source $SMD_ROOT/arp_scripts/setup_hosts_openmpi.sh
    fi
fi

export SMD_DEBUG_SMALLDATA_TIMING=1  # in psana/smalldata.py
export SMD_DEBUG_TIMING=1
export SMD_DEBUG_EVT_TIMING=1
export SMD_DEBUG_TIMING_EVERY=1
export SMD_DEBUG_DET=1
export SMD_DEBUG_DETOBJ=0
export SMD_DEBUG_AZAV_SHM=0
echo "[DEBUG] SMD_DEBUG_SMALLDATA_TIMING: $SMD_DEBUG_SMALLDATA_TIMING"
echo "[DEBUG] SMD_DEBUG_TIMING: $SMD_DEBUG_TIMING"
echo "[DEBUG] SMD_DEBUG_EVT_TIMING: $SMD_DEBUG_EVT_TIMING"
echo "[DEBUG] SMD_DEBUG_TIMING_EVERY: $SMD_DEBUG_TIMING_EVERY"
echo "[DEBUG] SMD_DEBUG_DET: $SMD_DEBUG_DET"
echo "[DEBUG] SMD_DEBUG_DETOBJ: $SMD_DEBUG_DETOBJ"
echo "[DEBUG] SMD_DEBUG_AZAV_SHM: $SMD_DEBUG_AZAV_SHM"
echo "Producer command: $SMD_ROOT/lcls2_producers/$PYTHONEXE $@"
if [ -n "$MPI_CORES" ]; then
    mpirun -n "$MPI_CORES" python -u -m mpi4py.run "$SMD_ROOT/lcls2_producers/$PYTHONEXE" "$@"
else
    mpirun python -u -m mpi4py.run "$SMD_ROOT/lcls2_producers/$PYTHONEXE" "$@"
fi
