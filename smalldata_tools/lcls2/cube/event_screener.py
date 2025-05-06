"""
Event screening and filtering system for LCLS2 data analysis.
This module provides a flexible framework for implementing various event filtering criteria
in LCLS2 data analysis workflows. It includes both simple filters for threshold, range,
and boolean checks, as well as composite filters for combining multiple filtering conditions.

Classes
-------
EventScreener : ABC
    Abstract base class defining the interface for all event screeners.
ThresholdFilter : EventScreener
    Filter events based on whether a value is above/below a threshold.
RangeFilter : EventScreener
    Filter events based on whether a value falls within a specified range.
BoolFilter : EventScreener
    Filter events based on boolean state matching.
CompositeFilter : EventScreener
    Combine multiple filters using logical AND/OR operations.

Examples 1
----------
>>> # Create a simple threshold filter
>>> threshold_filter = ThresholdFilter(detector, 'peak_value', threshold=100)
>>> passes, label = threshold_filter.apply(event)

Examples 2
----------
>>> # Create a composite filter with multiple conditions
>>> filters = [
...     ThresholdFilter(det1, 'value', threshold=10),
...     RangeFilter(det2, 'temp', min_value=20, max_value=30)
... ]
>>> composite = CompositeFilter(filters, require_all=True, label='valid_events')
>>> passes, label = composite.apply(event)

- All filters implement a common interface through the EventScreener ABC
- Filters can be labeled for tracking which conditions events satisfy
- Composite filters can combine multiple conditions using AND/OR logic
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Any, Optional, Union, Callable, Set, Tuple

import psana

from . import utils

logger = logging.getLogger(__name__)


class EventScreener(ABC):
    """
    Abstract base class for event screening/filtering in LCLS2 data analysis.
    The EventScreener provides a common interface for implementing various event filtering
    criteria. Each concrete implementation should define specific logic for determining
    whether an event passes certain conditions.

    Examples
    --------
    >>> class ThresholdScreener(EventScreener):
    ...     def passes(self, evt):
    ...         return self.get_scalar_value(evt) > threshold
    >>> screener = ThresholdScreener(detector, 'value', label='above_threshold')
    >>> passes, label = screener.apply(event)

    Notes
    -----
    Concrete implementations must override the `passes` method to define specific
    filtering criteria. The `apply` method provides a consistent interface for
    filter application across all implementations.
    """

    def __init__(self, detector: Any, data_key: str, label: Optional[str] = None):
        """
        Parameters
        ----------
        detector : Any
            Detector object that provides access to event data
        data_key : str
            Name of the data field to be accessed from the detector data
        label : Optional[str], default=None
            Optional label to identify this screener's results
        """
        self.detector = detector
        self.data_key = data_key
        self.label = label

    @abstractmethod
    def passes(self, value: Union[float, int, bool]) -> bool:
        """
        Check if the event passes this filter.

        Parameters
        ----------
        value : Union[float, int, bool]
            The value to be checked against the filter criteria

        Returns
        -------
        bool
            True if the event passes the filter, False otherwise
        """
        pass

    def apply(self, evt: psana.event.Event) -> Tuple[bool, Optional[str]]:
        """
        Apply the filter and return result and label if passed.

        Parameters
        ----------
        evt : psana.event.Event
            The event to be screened

        Returns
        -------
        Tuple[bool, Optional[str]]
            (passes_filter, label)
            If filter passes and has a label, returns (True, label)
            If filter passes but has no label, returns (True, None)
            If filter doesn't pass, returns (False, None)
        """
        value = self.get_scalar_value(evt)
        result = self.passes(value)
        if result and self.label:
            return True, self.label
        return result, None

    def get_scalar_value(self, evt: psana.event.Event) -> Union[float, int, bool]:
        """
        Return the scalar value from the detector data. It is assumed that the
        detector data is either a scalar or a dictionary with the data field as a key.

        If the detector does not return a scalar natively, this method must be
        overridden in the subclass to provide the appropriate logic.

        Parameters
        ----------
        evt : psana.event.Event
        """
        data = self.detector.data(evt)
        return data[self.data_key] if isinstance(data, dict) else data


class ThresholdFilter(EventScreener):
    def __init__(
        self,
        detector: Any,
        data_key: str,
        threshold: float,
        greater_than: bool = True,
        label: Optional[str] = None,
    ):
        """
        Filter that checks if a detector's value is above/below a threshold.

        Parameters
        ----------
        detector : Any
            Detector object that provides data to compare against threshold
        data_key : str
            Name of the data field to be accessed from the detector data
        threshold : float
            The threshold value to compare detector data against
        greater_than : bool, optional
            If True, passes events where detector value > threshold
            If False, passes events where detector value < threshold
            Default is True
        label : str, optional
            Optional label for the cube.
            Default is None
        """
        super().__init__(detector_name, data_key, label=label)
        self.threshold = threshold
        self.greater_than = greater_than

    def passes(self, value) -> bool:
        if value is None:
            return False

        if self.greater_than:
            return value > self.threshold
        else:
            return value < self.threshold

    def __str__(self) -> str:
        operator = ">" if self.greater_than else "<"
        label_str = f" [{self.label}]" if self.label else ""
        return f"{self.detector_name} {operator} {self.threshold}{label_str}"


class RangeFilter(EventScreener):
    """Filter that checks if a detector's value is within a specified range."""

    def __init__(
        self,
        detector: Any,
        data_key: str,
        min_value: float,
        max_value: float,
        label: Optional[str] = None,
    ):
        """
        Parameters
        ----------
        detector : Any
            Detector object that provides data to compare against range
        data_key : str
            Name of the data field to be accessed from the detector data
        min_value : float
            Minimum value of the range
        max_value : float
            Maximum value of the range
        """
        super().__init__(detector, data_key, label=label)
        self.min_value = min_value
        self.max_value = max_value

    def passes(self, value) -> bool:
        if value is None:
            return False
        return self.min_value <= value <= self.max_value

    def __str__(self) -> str:
        brackets = "[]"
        label_str = f" [{self.label}]" if self.label else ""
        return f"{self.detector.name} in {brackets[0]}{self.min_value}, {self.max_value}{brackets[1]}{label_str}"


class BoolFilter(EventScreener):
    """Filter that checks if a detector's boolean value matches the expected state."""

    def __init__(
        self,
        detector,
        data_key: str,
        expected_state: bool = True,
        label: Optional[str] = None,
    ):
        """
        Parameters
        ----------
        detector : Any
            Detector object that provides data to compare against expected state
        data_key : str
            Name of the data field to be accessed from the detector data
        expected_state : bool
            The expected boolean state to check against
        label : str, optional
            Optional label for the cube.
            Default is None
        """
        super().__init__(detector, data_key, label=label)
        self.expected_state = expected_state

    def passes(self, value) -> bool:
        if value is None:
            return False

        # Convert value to boolean if it's not already
        if not isinstance(value, bool):
            # For scalar values, treat non-zero as True
            try:
                value = bool(value)
            except (ValueError, TypeError):
                return False

        return value == self.expected_state

    def __str__(self) -> str:
        state_str = "True" if self.expected_state else "False"
        label_str = f" [{self.label}]" if self.label else ""
        return f"{self.detector.name} is {state_str}{label_str}"


class CompositeFilter(EventScreener):
    """
    Filter that combines multiple filters using logical operations.
    """

    def __init__(
        self,
        filters: List[EventScreener],
        require_all: bool = True,  # True: AND operation, False: OR operation
        label: Optional[str] = None,
    ):
        """
        Initialize a composite event screener.

        Parameters
        ----------
        filters : List[EventScreener]
            List of event screeners to be combined
        require_all : bool, optional
            If True, all filters must pass (AND operation).
            If False, any filter passing is sufficient (OR operation).
            Defaults to True.
        label : str, optional
            Label for this screener. Defaults to None. If None, it will concatenate the
            labels of the underlying filters that passes.
        """
        # Use a dummy key for the base class
        super().__init__(None, "composite", label=label)
        self.filters = filters
        self.require_all = require_all

    def get_scalar_value(self, evt: psana.event.Event) -> None:
        """Override to prevent using the base class implementation."""
        return None

    def passes(self, value: Any) -> bool:
        """
        This is a placeholder to satisfy the abstract method requirement.
        The actual logic is in apply().
        """
        return False  # Never called directly

    def apply(self, evt: psana.event.Event) -> Tuple[bool, Optional[str]]:
        """
        Apply the filter and determine the appropriate label.

        Parameters
        ----------
        evt : psana.event.Event
            The event to filter

        Returns
        -------
        Tuple[bool, Optional[str]]
            (passes_filter, label)
        """
        # Check if composite filter passes
        passing_labels = []
        if self.require_all:
            # AND logic: all filters must pass
            for f in self.filters:
                passes, label = f.apply(evt)
                logger.debug(f"Filter: {f}, passes: {passes}, label: {label}")
                if not passes:
                    return False, None
                if label:
                    passing_labels.append(label)
        else:
            # OR logic: at least one filter must pass
            any_pass = False
            for f in self.filters:
                passes, label = f.apply(evt)
                if passes:
                    any_pass = True
                    if label:
                        passing_labels.append(label)

            if not any_pass:
                return False, None

        # If this composite has a label, use it
        if self.label:
            return True, self.label

        # Otherwise, concatenate labels from passing filters
        if passing_labels:
            return True, "_".join(passing_labels)
        else:
            return True, None

    def __str__(self) -> str:
        """Return a string representation of the composite filter."""
        operation = "AND" if self.require_all else "OR"
        filter_strs = [str(f) for f in self.filters]
        label_str = f" [{self.label}]" if self.label else ""
        return f"({' {operation} '.join(filter_strs)}){label_str}"
