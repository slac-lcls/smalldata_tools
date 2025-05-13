from enum import Enum  # , StrEnum py3.11)
import logging
import time


logger = logging.getLogger(__name__)


class PsanaNode(str, Enum):
    SMD0 = "SMD0"
    EB = "EB"
    BD = "BD"
    SRV = "SRV"


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
