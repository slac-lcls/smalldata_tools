import sys
import os
import logging
import time
import psana
from enum import Enum  # , StrEnum py3.11)
from pathlib import Path
from mpi4py import MPI

COMM = MPI.COMM_WORLD
rank = COMM.Get_rank()
size = COMM.Get_size()


logger = logging.getLogger(__name__)


class PsanaNode(str, Enum):
    SMD0 = "SMD0"
    EB = "EB"
    BD = "BD"
    SRV = "SRV"


def get_psana_node_type(ds: psana.DataSource) -> PsanaNode:
    if rank == 0:
        return PsanaNode.SMD0

    if ds.is_srv():
        print(f"SRV World Rank: {rank}")
        psana_node = PsanaNode.SRV
        psana_comm = MPI.COMM_NULL
        psana_rank = MPI.UNDEFINED
        bd_rank = MPI.UNDEFINED
        return PsanaNode.SRV

    comms = ds.comms

    psana_comm = comms.psana_comm
    psana_size = psana_comm.Get_size()
    psana_rank = psana_comm.Get_rank()

    bd_comm = comms.bd_comm
    if bd_comm is not None:
        bd_size = bd_comm.Get_size()
        bd_rank = bd_comm.Get_rank()
    else:
        bd_rank = MPI.UNDEFINED

    if bd_rank == 0:
        psana_node = PsanaNode.EB
    elif bd_rank > 0:
        psana_node = PsanaNode.BD

    if ds.unique_user_rank():
        print(f"MPI COMM sizes:\tWorld: {size},\tpsana: {psana_size},\tbd: {bd_size}\n")
        n_eb = int(os.environ.get("PS_EB_NODES"))
        n_srv = int(os.environ.get("PS_SRV_NODES"))
        n_bd = size - n_eb - n_srv - 1
        print(
            f"PSANA nodes: # EB nodes (total): {n_eb}, # BD nodes (total): {n_bd}, "
            f"# BD nodes per EB: {n_bd / n_eb}\n\n"
        )
    return psana_node


def get_config_file(name, folder_path):
    """Find config file using pathlib"""
    folder = Path(folder_path)
    target_file = folder / f"cube_config_{name}.py"

    if target_file.exists():
        return target_file.stem  # return the file name without extension
    else:
        if rank == 0:
            logger.error(f"Config file {target_file} not found.")
        sys.exit(1)


def reduce_by_key(reduced: dict, new_data: dict, reduction_func: callable) -> dict:
    """
    Recursively reduce the new data by the keys of the reduced data using the provided
    reduction function.

    Parameters
    ----------
    reduced : dict
        The dictionary holding the reduced data.
    new_data : dict
        The new data to be added to the reduced dictionary.
    reduction_func : callable
        The function to reduce the data. It should take two arguments: the existing
        value and the new value.

    Returns
    -------
    dict
        The reduced dictionary.
    """
    for key, value in new_data.items():
        if isinstance(value, dict):
            if key not in reduced:
                reduced[key] = {}
            reduce_by_key(reduced[key], new_data[key], reduction_func)
        else:
            if key in reduced:
                reduced[key] = reduction_func(reduced[key], value)
            else:
                reduced[key] = value
    return reduced


def append_reduction(reduced, new):
    if isinstance(reduced, list):
        reduced.append(new)
    else:
        reduced = [reduced]
        reduced.append(new)
    return reduced


def timing_decorator(func):
    """
    Decorator to measure execution time of methods.

    Example:
    --------
        @timing_decorator
        def apply(self, evt):
            ...
    """

    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        logger.debug(f"{func.__name__} took {(end-start)*1000:.2f} ms")
        return result

    return wrapper
