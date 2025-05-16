import os
import numpy as np
import h5py
import logging

from typing import Union, List, Any, Optional

logger = logging.getLogger(__name__)


def append_to_dataset(h5_dataset, data, is_sum_dataset=False):
    """
    Append data to an existing HDF5 dataset or sum data to it.

    Parameters
    ----------
    h5_dataset : h5py.Dataset
        The HDF5 dataset to which the data will be appended.
        Must be resizable (created with maxshape=(None, ...)) if appending.
    data : array-like
        The data to append to the dataset.
    is_sum_dataset : bool, optional
        If True, the data will be summed instead of appended with the existing dataset.
        Default is False.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the shape of data is incompatible with the dataset.
    RuntimeError
        If the dataset is not resizable when trying to append.
    """
    # Convert to numpy array if not already
    data_array = np.asarray(data)

    logger.debug(
        f"###  Add to dataset: {h5_dataset.name}, shape: {h5_dataset.shape}, "
        f"maxshape: {h5_dataset.maxshape}, "
        f"data: shape: {data_array.shape}"
    )

    if is_sum_dataset:
        # For datasets marked for summing
        try:
            current_value = h5_dataset[
                ()
            ]  # Fixed variable name from dataset to h5_dataset
            h5_dataset[()] = current_value + data_array
        except ValueError as e:
            # Handle shape mismatch
            raise ValueError(
                f"Cannot sum data with shape {data_array.shape} to dataset "
                f"with shape {h5_dataset.shape}: {e}"
            )
    else:
        # Append data to the dataset
        if data_array.ndim == 0:
            # If data is a scalar, convert to 1D array
            data_array = np.expand_dims(data_array, 0)

        # Check if dataset is resizable
        if (
            h5_dataset.maxshape[0] is None
            or h5_dataset.maxshape[0] > h5_dataset.shape[0]
        ):
            # Check if shapes match in all dimensions except the first one
            if h5_dataset.ndim > 1:
                if h5_dataset.shape[1:] != data_array.shape:
                    raise ValueError(
                        f"Shape mismatch: dataset {h5_dataset.name} dimensions "
                        f"{h5_dataset.shape[1:]} vs data dimensions {data_array.shape}"
                    )
            else:
                # For 1D datasets, we just append without shape checks
                pass

            # Add a dimension to the data for appending
            data_array = np.expand_dims(data_array, 0)

            # Resize only the first dimension
            current_size = h5_dataset.shape[0]
            new_size = current_size + len(data_array)
            h5_dataset.resize((new_size,) + h5_dataset.shape[1:])

            # Add new data at the end of the first dimension
            h5_dataset[current_size:new_size] = data_array
        else:
            raise RuntimeError(
                f"Dataset {h5_dataset.name} is not resizable along the first dimension"
            )


def create_dataset_from_key_value(
    parent: Union[h5py.File, h5py.Group],
    key: str,
    value: Any,
    sum_dataset: List[str] = [],
    len_data: Optional[int] = None,
    add_at_index: Optional[int] = None,
) -> h5py.Dataset:
    """
    Create a new dataset in the HDF5 file or group.

    Parameters
    ----------
    parent : Union[h5py.File, h5py.Group]
        The HDF5 File or Group object where the dataset will be created.
    key : str
        The name/key of the dataset.
    value : Any
        The value to store in the dataset.
    sum_dataset : List[str], optional
        List of keys for which to create scalar datasets for summing. Default is ['count'].
    len_data : Optional[int], optional
        Initial size of the first dimension for the dataset.
    add_at_index : Optional[int], optional
        Index at which to add data in the dataset. If None, use data to create data set.

    Returns
    -------
    h5py.Dataset
        The newly created dataset.
    """
    if add_at_index is not None and len_data is None:
        raise ValueError("len_data must be specified if add_at_index is used")

    logger.info(f"### Create dataset in {parent} with name: {key}")

    # Create new dataset with the dict data
    data = np.asarray(value)

    if key in sum_dataset:  # New data will be summed
        parent.create_dataset(key, data=data)

    else:  # New data will be appended
        """Note: if data in dict are multiple entries, then the first index is
        the "dataset" index. But this is not the case for the cube, bins are
        passed one at a time. So let's add a dimension to the data so that we
        can append to it when new data for this dataset comes."""
        if data.ndim == 0:
            # If data is a scalar, convert to 1D array
            data = np.expand_dims(data, 0)
        data = np.expand_dims(data, 0)  # Add a new dimension for appending
        dataset = parent.create_dataset(
            key, data=data, maxshape=(None,) + data.shape[1:]
        )

        if False:
            # OTHER APPROACH. Can be modified for add_at_index

            # Create empty dataset with shape info from the array in data_dict
            data = np.asarray(value)

            # Get the shape and dtype from the sample data
            sample_shape = data.shape
            dtype = data.dtype

            # Create an empty, resizable dataset
            if len(sample_shape) > 0:
                # Multi-dimensional array (potentially)
                if len_data is not None:
                    initial_shape = (len_data,) + sample_shape[1:]
                else:
                    initial_shape = (0,) + sample_shape[
                        1:
                    ]  # Start with 0 in first dimension
                max_shape = (None,) + sample_shape[1:]  # Unlimited first dimension
            else:
                # Scalar or 1D array
                if len_data is not None:
                    initial_shape = (len_data,)
                else:
                    initial_shape = (0,)
                max_shape = (None,)

            dataset = parent.create_dataset(
                key,
                shape=initial_shape,
                maxshape=max_shape,
                dtype=dtype,
                chunks=True,  # Let h5py choose an appropriate chunk size
            )

        return dataset


def add_dict_to_h5(
    parent: Union[h5py.File, h5py.Group],
    data_dict: dict,
    sum_dataset: List[str] = [],
    len_data: Optional[int] = None,
    add_at_index: Optional[int] = None,
) -> None:
    """
    Add data to an HDF5 file or group, creating groups and datasets as needed. For now this
    function is intented for datasets that will be appended one-by-one. Future versions will
    allow for adding data at a specific index and summing datasets.

    Parameters
    ----------
    parent : Union[h5py.File, h5py.Group]
        The HDF5 File or Group object to write to.
    data_dict : dict
        Nested dictionary containing the data structure to write.
    sum_dataset : List[str], optional
        Keys for which to sum instead of appending data. Default is ['count'].
    len_data : Optional[int], optional
        Size of the first dimension for all datasets.
        Comment: Could use the add_at_index as the potential size of the dataset, and
                 resize if needed...
    add_at_index : Optional[int], optional
        Index at which to add data in datasets. If None, appends data.
        len_data must be specified if add_at_index is used.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If add_at_index is specified without len_data
    """
    if add_at_index is not None and len_data is None:
        raise ValueError("len_data must be specified if add_at_index is used")

    for key, value in data_dict.items():
        if isinstance(value, dict):
            # If group exists, use it; otherwise create it
            if key in parent:
                group = parent[key]
            else:
                group = parent.create_group(key)

            # Always recurse into the group
            add_dict_to_h5(
                group,
                value,
                sum_dataset=sum_dataset,
                len_data=len_data,
                add_at_index=add_at_index,
            )
        else:
            # If dataset already exists, append data or sum, else create it
            if key in parent:
                if add_at_index is not None:
                    # TODO: Make function to add data at specified index
                    raise NotImplementedError(
                        "Adding data at a specific index is not \
                                               implemented yet."
                    )
                else:
                    # Append data to existing dataset
                    append_to_dataset(parent[key], value)

            else:
                dataset = create_dataset_from_key_value(
                    parent,
                    key,
                    value,
                    sum_dataset=sum_dataset,
                    len_data=len_data,
                    add_at_index=add_at_index,
                )


if __name__ == "__main__":
    """
    THIS SHOULD PROBABLY BE TURNED INTO PYTEST TESTS

    Example script to demonstrate the usage of h5_utils module for managing HDF5 file
    operations.
    This script shows how to:
        - Create nested data structures using dictionaries with numpy arrays
        - Write data to HDF5 files using add_bin_dict_to_h5()
        - Append new data to existing datasets
        - Handle different data shapes and types
        - Print dataset information

    Example data structures are created with:
        - Scalar values
        - 1D arrays
        - 2D arrays
        - Nested dictionaries

    This serves as a practical example of the core functionality provided by h5_utils
    module.
    """

    def print_dataset_info(ds):
        print(f"{ds.name}:")
        print(f"shape: {ds.shape}")
        print(f"maxshape: {ds.maxshape}")
        print(f"dtype: {ds.dtype}")
        print(f"data:\n{ds[()]}")
        print("\n")

    # Build some example data
    bin_data_0 = {
        "d1": {"s1": np.arange(5), "s2": np.ones((3, 6))},
        "d2": {"s1": 10, "s2": np.arange(10, 20)},
    }

    bin_data_1 = {
        "d1": {"s1": 2 * np.arange(5), "s2": 2 * np.ones((3, 6))},
        "d2": {"s1": 11, "s2": np.arange(30, 40)},
        "d3": np.random.randint(10, size=(6, 11)),
    }

    f = h5py.File("test_h5_utils.h5", "w")
    add_bin_dict_to_h5(f, bin_data_0)

    print("\n")
    print_dataset_info(f["d1"]["s1"])
    print_dataset_info(f["d1"]["s2"])
    print_dataset_info(f["d2"]["s1"])
    print_dataset_info(f["d2"]["s2"])

    print("\nAdd new data")
    add_bin_dict_to_h5(f, bin_data_1)

    print_dataset_info(f["d1"]["s1"])
    print_dataset_info(f["d1"]["s2"])
    print_dataset_info(f["d2"]["s1"])
    print_dataset_info(f["d2"]["s2"])
    print_dataset_info(f["d3"])

    f.close()
    os.remove(f)
