#!/bin/bash

# Export useful path
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
SMD_ROOT=`echo $MYDIR | sed  "s|/arp_scripts||g"`

PY_EXE="$SMD_ROOT/smalldata_tools/ana_funcs/svd_waveform/make_waveform_basis_psana2.py"

echo "Sourcing LCLS-II environment"
export SIT_ENV_DIR="/sdf/group/lcls/ds/ana"
source $SIT_ENV_DIR/sw/conda2/manage/bin/psconda.sh
PYTHONPATH="$SMD_ROOT:$PYTHONPATH"

python $PY_EXE $*