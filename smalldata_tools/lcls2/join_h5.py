#!/sdf/group/lcls/ds/ana/sw/conda2/inst/envs/ps-4.6.3/bin/python

# Courstesy of cpo

import sys
import h5py
import glob
from psana.smalldata import _get_missing_value, is_unaligned


def myjoin(joined_file):
    # Get part files list
    # Assumes that the files are of the form <joined_file)_part<N>.h5
    joined_file = str(joined_file)
    joined_file = h5py.File(joined_filename, "w", libver="latest")
    files = glob.glob(joined_filename.replace(".h5", "_part*.h5"))

    # discover all the dataset names
    file_dsets = {}

    def assign_dset_info(name, obj):
        # TODO check if name contains unaligned, if so ignore
        if isinstance(obj, h5py.Dataset):
            tmp_dsets[obj.name] = (obj.dtype, obj.shape)

    all_dsets = []
    for fn in files:
        tmp_dsets = {}
        f = h5py.File(fn, "r")
        f.visititems(assign_dset_info)
        file_dsets[fn] = tmp_dsets
        all_dsets += list(tmp_dsets.keys())
        f.close()

    all_dsets = set(all_dsets)

    # h5py requires you declare the size of the VDS at creation
    # (we have been told by Quincey Koziol that this is not
    # necessary for the C++ version).
    # so: we must first loop over each file to find # events
    #     that come from each file
    # then: do a second loop to join the data together

    for dset_name in all_dsets:

        # part (1) : loop over all files and get the total number
        # of events for this dataset
        total_events = 0
        for fn in files:
            dsets = file_dsets[fn]
            if dset_name in dsets.keys():
                dtype, shape = dsets[dset_name]
                total_events += shape[0]

            # this happens if a dataset is completely missing in a file.
            # to maintain alignment, we need to extend the length by the
            # appropriate number and it will be filled in with the
            # "fillvalue" argument below.  if it's unaligned, then
            # we don't need to extend it at all.
            elif not is_unaligned(dset_name):
                if "/timestamp" in dsets:
                    total_events += dsets["/timestamp"][1][0]

        combined_shape = (total_events,) + shape[1:]

        layout = h5py.VirtualLayout(shape=combined_shape, dtype=dtype)

        # part (2): now that the number of events is known for this
        # dataset, fill in the "soft link" that points from the
        # master file to all the smaller files.
        index_of_last_fill = 0
        for fn in files:

            dsets = file_dsets[fn]

            if dset_name in dsets.keys():
                _, shape = dsets[dset_name]
                vsource = h5py.VirtualSource(fn, dset_name, shape=shape)
                layout[index_of_last_fill : index_of_last_fill + shape[0], ...] = (
                    vsource
                )
                index_of_last_fill += shape[0]

            else:
                # only need to pad aligned data with "fillvalue" argument below
                if is_unaligned(dset_name):
                    pass
                else:
                    if "/timestamp" in dsets:
                        n_timestamps = dsets["/timestamp"][1][0]
                        index_of_last_fill += n_timestamps

        joined_file.create_virtual_dataset(
            dset_name, layout, fillvalue=_get_missing_value(dtype)
        )

    joined_file.close()
    return 0
