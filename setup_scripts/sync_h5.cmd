#!/bin/bash
rsync -avu $FFB_BASE/hdf5/smalldata/*.h5 $PSANA_BASE/hdf5/smalldata
rsync -avu $FFB_BASE/hdf5/cube/*.h5 $PSANA_BASE/hdf5/cube
