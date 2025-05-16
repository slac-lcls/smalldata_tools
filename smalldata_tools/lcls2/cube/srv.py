import os
import logging
import h5py

import psana

from enum import Enum  # , StrEnum py3.11)
from typing import Union, Any
from dataclasses import dataclass, field
from pathlib import Path

from . import h5_utils
from . import utils

from mpi4py import MPI

COMM = MPI.COMM_WORLD
rank = COMM.Get_rank()
size = COMM.Get_size()

logger = logging.getLogger(__name__)


SIT_PSDM_DATA = Path(os.environ.get("SIT_PSDM_DATA", "/sdf/data/lcls/ds/"))


class SrvMsgType(Enum):
    """
    Message types for the server.
    """

    DONE = 0
    NEW_BIN = 1
    STEP_DONE = 2


@dataclass
class BinData:
    """
    This class stores the data of a worker for a given bin.
    It is also used on the SRV to store the aggregated data for a given bin.

    Attributes:
    ----------
        cube_label (str):
            Label identifying the cube this data belongs to.
            Ex: "laser_on"
        data (dict):
            Dictionary containing the bin data
        bin_index (int):
            Index identifying which bin this data belongs to
        bin_info (Any):
            Additional metadata or information about the bin.
            Can be of any type depending on the binning method used.

    """

    cube_label: str
    bin_index: int
    data: dict
    bin_info: dict = None

    def __add__(self, other):
        """
        Add two BinData objects together.
        This method will sum the data dictionaries of the two objects.

        Note: this method will modify the original self.data dictionary in place.
        """
        if self.cube_label != other.cube_label:
            raise ValueError("Cannot add BinData with different cube labels.")

        if self.bin_index != other.bin_index:
            raise ValueError("Cannot add BinData with different bin indices.")

        if self.bin_info != other.bin_info:
            raise ValueError("Cannot add BinData with different bin info.")

        new_data = utils.reduce_by_key(self.data, other.data, lambda x, y: x + y)
        return BinData(
            self.cube_label, self.bin_index, new_data, bin_info=self.bin_info
        )


@dataclass
class SrvCubeMessage:
    """
    Message to be sent to the server.
    """

    msg_type: SrvMsgType
    sender: int
    payload: Union[dict, BinData] = field(default_factory=dict)


class CubeSrv:
    """
    Class for handling data collection from MPI workers and writing to HDF5.
    This class receives binned data from worker nodes (BD nodes) and combines partial
    results before writing the final data to an HDF5 file.

    The class operates in two modes:
    - Scan mode: Each BD worker processes all bins
    - Non-scan mode: Each bin is processed by a single BD worker (TODO: NOT IMPLEMENTED)

    Attributes:
        n_bd (int): Number of BD (Bin Data) nodes that will send data
        file_handle (h5py.File): Handle to the HDF5 file for writing data

    Methods:
        set_file_handle: Sets up the HDF5 file for writing data
        recv: Receives messages from BD nodes
        yield_from_bd: Yields bin data from BD nodes until all workers are done
        run: Main loop that collects and processes data from workers
        process_bin: Processes a complete bin and writes to HDF5
        write_bin_to_h5: Writes bin data to the HDF5 file

    Args:
        file_kwargs (dict): Keywords arguments for file handling
        scan_mode (bool, optional): If True, each BD processes all bins. Defaults to True.
    """

    def __init__(self, scan_mode: bool = True):
        n_eb = int(os.environ.get("PS_EB_NODES"))
        n_srv = int(os.environ.get("PS_SRV_NODES"))

        # Number of BD nodes working on a given bin.
        # In scan mode, a bin is distributed on all BD. Else 1 BD/bin:
        self.n_bd = size - n_eb - n_srv - 1 if scan_mode else 1
        logger.info(f"Server expecting data from {self.n_bd} BD nodes.")
        self.file_handle = None

    def set_file_handle(self, run_num: int, exp: str, filepath: str = None):
        """
        Get the file handle to write data to.
        """
        if filepath is not None:
            logger.info(f"User custom filepath {filepath}.")
            filepath = Path(filepath)
        else:
            filepath = Path(SIT_PSDM_DATA) / exp[:3] / exp / "hdf5/smalldata/cube"
            logger.info(f"Using default filepath {filepath}.")

        filename = f"cube_{exp}_r{run_num:04d}.h5"
        filename = filepath / filename
        logger.info(f"Write cube data to {filename}.")
        self.file_handle = h5py.File(filename, "w")

    def recv(self):
        """
        Get binned data from the BD nodes.
        This method will block until a message is received.
        """
        status = MPI.Status()
        msg = COMM.recv(source=MPI.ANY_SOURCE, status=status)
        sender = status.Get_source()
        logger.debug(f"Received msg from {sender}: {msg.msg_type}\n")
        return sender, msg

    def yield_from_bd(self):
        all_done = False
        count_done = 0
        while not all_done:
            sender, msg = self.recv()
            if msg.msg_type == SrvMsgType.DONE:
                count_done += 1
                if count_done == size - 1:
                    all_done = True
                    logger.info("All psana nodes are done.")
            else:
                yield msg

    def run(self):
        # Dictionary to store partially combined data
        bin_cache: Dict[Tuple[str, int], BinData] = {}
        # keep track of how many BD contributed to a bin:
        bin_count: Dict[Tuple[str, int], int] = {}
        # keep track of bd who finished a step:
        step_done: Dict[int, int] = {}

        # for bin_data in self.yield_from_bd():
        for msg in self.yield_from_bd():
            if msg.msg_type == SrvMsgType.NEW_BIN:
                bin_data = msg.payload
                key = (bin_data.cube_label, bin_data.bin_index)
                logger.debug(f"Got bin data for {key}")

                if key in bin_cache:
                    # Combine with existing data
                    bin_cache[key] = bin_cache[key] + bin_data
                    bin_count[key] += 1
                else:
                    # Store new data
                    bin_cache[key] = bin_data
                    bin_count[key] = 1
            elif msg.msg_type == SrvMsgType.STEP_DONE:
                step_idx = msg.payload["step_value"] - 1  # step starts from 1
                if step_idx in step_done:
                    step_done[step_idx] += 1
                else:
                    # First time we see this step
                    step_done[step_idx] = 1

            to_be_deleted = []
            for key in bin_cache.keys():
                logger.debug(f"key: {key}, count: {bin_count[key]}")
                if key[1] in step_done:
                    if step_done[key[1]] == self.n_bd:
                        """
                        Note: we cannot rely on the bin_count here, because in some cases, not all
                        BD will send data for a given bin. This is the case for example when there
                        are more bd than there are batches of events in a step. In this case, the
                        bin_count will be less than n_bd. Relying on the step_done count is safer.
                        """
                        # All data for this bin has been received
                        logger.info(
                            f"All data for {key} received. Number of BD who contributed: {bin_count[key]}"
                        )
                        self.process_bin(bin_cache[key])
                        to_be_deleted.append(key)
            for key in to_be_deleted:
                del bin_cache[key]
                # del bin_count[key]  # Let's keep the count for debugging

        logger.info("Server done, close file handle.")
        self.file_handle.close()

    def process_bin(self, bin_data: BinData):
        logger.info(
            f"Processing bin for cube label: {bin_data.cube_label}, bin index: {bin_data.bin_index}"
        )
        logger.info(f"Bin info: {bin_data.bin_info}")
        # TODO: Add processing pipeline on the binned data here.
        step_info = bin_data.bin_info
        # TODO: add docstring to h5. For now it can't handle strings
        step_info.pop("step_docstring")
        h5_utils.add_dict_to_h5(
            parent=self.file_handle,
            data_dict=step_info,
        )
        self.write_bin_to_h5(bin_data)

    def write_bin_to_h5(self, bin_data: BinData):
        """
        Write the bin data to the HDF5 file.
        """
        if bin_data.cube_label is not None:
            # Make cube group if first time we see this label
            if bin_data.cube_label not in self.file_handle:
                self.file_handle.create_group(bin_data.cube_label)
            handle = self.file_handle[bin_data.cube_label]
        else:
            handle = self.file_handle

        logger.debug(f"Use handle {handle} to write bin data")
        h5_utils.add_dict_to_h5(
            parent=handle,
            data_dict=bin_data.data,
        )
