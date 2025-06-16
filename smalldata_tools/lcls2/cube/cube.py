from abc import ABCMeta, abstractmethod
import logging
import psana

from typing import Any, Callable
from mpi4py import MPI

COMM = MPI.COMM_WORLD
rank = COMM.Get_rank()
size = COMM.Get_size()

from . import event_engine
from .utils import PsanaNode, get_psana_node_type, reduce_by_key, append_reduction
from .srv import SrvMsgType, SrvCubeMessage, BinData
from . import event_screener

logger = logging.getLogger(__name__)


class Cube(metaclass=ABCMeta):
    def __init__(
        self,
        run: psana.psexp.run.Run,
        engine: Callable = event_engine.smalldata_tools_engine,
        processors: list = None,
        event_screener: event_screener.EventScreener = None,
        **kwargs: dict,
    ):
        """
        Initialize a Cube instance for event processing.
        This class handles event processing for LCLS2 data using a specified processing engine
        and a list of event processors to be passed to the engine at each event. Typically each processor handles
        a single detector.
        The processing engine is responsible for handling the event data and applying the specified
        processors to the event data. By default, we assume that smalldata_tools is used to process detectors data.

        Parameters
        ----------
        run : psana.psexp.run.Run
            psana run object containing the experiment data to be processed.
        engine : callable, optional
            The processing engine to be used for event handling. Defaults to smalldata_tools_engine.
        processors : list, optional
            List of processor objects to be applied to each event. Defaults to empty list.
        Additional keyword arguments:
            bin_processors : list, optional
                List of processor objects to be applied to each bin after processing all events.

        """
        self._run = run
        self.engine = engine
        self.processors = processors if processors is not None else []
        self.event_screener = event_screener
        self.bin_processors = kwargs.get("bin_processors", [])
        # Reduction passed to utils.reduce_by_key to combine a bin from different
        # BD (simple sum):
        self.reduction_func = lambda x, y: x + y
        self.node_type = get_psana_node_type(run.ds)

    def add_processors(self, processors):
        """
        Add one or more processors to the cube's list of processors to run on each event.
        This method extends the existing list of processors by adding new processor(s).

        Parameters
        ----------
        processors : object or list
            A single processor object or a list of processor objects to be added

        Returns
        -------
        None

        Examples
        --------
        >>> cube.add_processors(my_processor)
        >>> cube.add_processors([processor1, processor2])
        """

        if not isinstance(processors, list):
            processors = [processors]
        self.processors.extend(processors)

    def set_event_screener(self, screener: event_screener.EventScreener):
        self.event_screener = screener

    def add_bin_processors(self, processors):
        """
        Add one or more processors to the cube's list of processors at the end of each bin.
        This method extends the existing list of processors by adding new processor(s).

        Parameters
        ----------
        processors : object or list
            A single processor object or a list of processor objects to be added

        Returns
        -------
        None

        Examples
        --------
        >>> cube.add_processors(my_processor)
        >>> cube.add_processors([processor1, processor2])
        """

        if not isinstance(processors, list):
            processors = [processors]
        self.bin_processors.extend(processors)

    def process_event(self, evt):
        evt_data = {}
        for proc in self.processors:
            proc_data = {}
            proc_data = self.engine(evt, proc)
            evt_data.update(proc_data)
        return evt_data

    def reduce_bin_data(self, binned_data, new_data):
        """
        Reduce cube data by applying a reduction function recursively based on binned data structure.

        Parameters
        ----------
        data : dict
            The input data to be reduced according to the cube's binning structure.

        Returns
        -------
        array-like
            The reduced data following the cube's binning structure and reduction function.

        Notes
        -----
        This method uses the cube's pre-defined binned_data structure and reduction_func
        to recursively reduce the input data. The reduction is performed by key at each
        level of the binned data hierarchy.
        """
        return reduce_by_key(binned_data, new_data, self.reduction_func)

    def send_bin(self, binned_data: dict, dest: int = size - 1) -> None:
        """
        Send the binned data to the server node.

        Parameters
        ----------
        binned_data : dict
            The binned data to be sent to the server node.

        Returns
        -------
        None

        Notes
        -----
        This method is intended to be called on the backend nodes (BD) to send the
        binned data to the server node (SRV) for further processing or storage.
        """
        if COMM != MPI.COMM_NULL:
            msg = SrvCubeMessage(
                msg_type=SrvMsgType.NEW_BIN, sender=rank, payload=binned_data
            )
            COMM.send(msg, dest=dest)  # The last rank is the srv node

    @abstractmethod
    def run(self):
        """
        Abstract method to be implemented by subclasses.
        This method should define how the event loop is run for the specific
        cube type (e.g., fly scan or step scan).
        """
        pass


class CubeFlyScan(Cube):
    def run(self):
        for nevt, evt in enumerate(self._run.events()):
            if nevt % 500 == 0:
                print(f"rank {rank}: {nevt}")
            evt_data = {}
            evt_data = self.process_event(evt)


class CubeStepScan(Cube):
    def __init__(
        self,
        run: psana.psexp.run.Run,
        engine: callable = event_engine.smalldata_tools_engine,
        processors: list = None,
        event_screener: event_screener.EventScreener = None,
        **kwargs: dict,
    ):
        super().__init__(run, engine, processors, event_screener, **kwargs)
        self.step = None
        self.step_data = {}

    def get_scan_detector(self) -> dict:
        """
        Get the step detector from the run object.
        This method is a placeholder and should be implemented in subclasses
        to define how the step detector is obtained for the specific cube type.
        """
        scanlist = [k[0] for k in self._run.scaninfo]
        scan_dets = {}
        for scan in scanlist:
            scan_dets[scan] = self._run.Detector(scan)
        return scan_dets

    def run(self) -> None:
        scan_dets = self.get_scan_detector()
        step_data = {}
        nstep = -1

        for nstep, step in enumerate(self._run.steps()):
            logger.debug(f"Rank {rank}: Begin step: {nstep}")

            # Get the step scan data
            step_data = {name: scan_det(step) for name, scan_det in scan_dets.items()}

            logger.debug(f"Step values: {step_data}")
            binned_data = {}

            # Loop over events in the step
            nevt = -1
            for nevt, evt in enumerate(step.events()):
                # Check if event passes filters and get potential cube label
                cube_label = None
                if self.event_screener:
                    passes, cube_label = self.event_screener.apply(evt)
                    if not passes:
                        continue

                # Get and process the event data
                evt_data = {}
                evt_data = self.process_event(evt)

                # Format the binned data according to the cube label
                # Note that if the label is None, the key is None. This is fine and can
                # be handled like any other key.
                if cube_label not in binned_data:
                    binned_data[cube_label] = {}
                binned_data[cube_label] = self.reduce_bin_data(
                    binned_data[cube_label], evt_data
                )

            # Format and send the binned data to the srv node. Each label is sent
            # separately.
            for label, binned_dict in binned_data.items():
                logger.debug(f"Formatting binned data for label {label}")
                bin_data = BinData(
                    cube_label=label,
                    bin_index=step_data["step_value"] - 1,  # step value starts at 1
                    bin_info=step_data,
                    data=binned_dict,
                )
                self.send_bin(bin_data, dest=size - 1)
            # Send the step done flag to the server node:
            logger.debug(f"Rank {rank}: number of events in step {nstep}: {nevt}")
            COMM.send(
                SrvCubeMessage(
                    msg_type=SrvMsgType.STEP_DONE, sender=rank, payload=step_data
                ),
                dest=size - 1,
            )

        if self.node_type == PsanaNode.SMD0:
            msg = SrvCubeMessage(msg_type=SrvMsgType.SMD0_DONE, sender=rank)
        elif self.node_type == PsanaNode.EB:
            msg = SrvCubeMessage(msg_type=SrvMsgType.EB_DONE, sender=rank)
        elif self.node_type == PsanaNode.BD:
            msg = SrvCubeMessage(
                msg_type=SrvMsgType.BD_DONE, sender=rank, payload=nstep
            )
        else:
            logger.error(f"Unexpected psana node type: {self.node_type}")
        COMM.send(msg, dest=size - 1)
        logger.info(f"Rank {rank} done at after {nstep} steps.")


def get_cube(
    run: psana.psexp.run.Run,
    engine: callable = event_engine.smalldata_tools_engine,
    **kwargs: dict,
) -> Cube:
    """
    Factory function to create a Cube object based on run type.
    This function determines whether the input run is a step scan or fly scan and
    returns the appropriate Cube object instance.

    Parameters
    ----------
    run : psana.psexp.run.Run
        The psana run object to create a cube from
    engine : callable, optional
        The engine to use for processing events, defaults to smalldata_tools_engine
    kwargs : dict, optional
        Additional keyword arguments to be passed to the Cube constructor. See the Cube
        class for details. This can include processors, event_screener, etc.

    Returns
    -------
    Cube
        Either a CubeFlyScan or CubeStepScan instance depending on the run type
    """
    if run.scaninfo == {}:
        is_scan = False
        if run.ds.unique_user_rank():
            logger.info("Instantiating CubeFlyScan.")
            logger.info("Not implemented yet.")
        raise NotImplementedError
        #cube = CubeFlyScan(run, engine=engine, **kwargs)
    else:
        is_scan = True
        if run.ds.unique_user_rank():
            logger.info("Instantiating CubeStepScan.")
        cube = CubeStepScan(run, engine=engine, **kwargs)
    return cube
