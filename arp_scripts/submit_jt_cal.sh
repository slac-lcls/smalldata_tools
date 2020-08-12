#!/bin/bash

source /reg/g/pcds/engineering_tools/xpp/scripts/pcds_conda
ABS_PATH=/reg/g/psdm/sw/tools/smalldata_tools/examples

sbatch --nodes=2 --time=5 $ABS_PATH/jt_cal.py "$@"
