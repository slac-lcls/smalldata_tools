"""
Collections of event engine tools for LCLS2. These engines defines how to process
a given event for a given processor. It is meant to more easily support various
event-base implementation of detectors processing in psana.

The engine must return data in a nested dictionary format, where the key is the processor
name or detector identifier, and the value is the corresponding data, in an int, float, array,
list or dictionary format.
If in a dictionary form, it is recommended to add a count key to keep track of
the number of events processed for that processor just before returning the data: 
data.update['count'] = 1
"""

import sys
import logging

logger = logging.getLogger(__name__)


try:
    from smalldata_tools.common.detector_base import getUserData, getUserEnvData

    def smalldata_tools_engine(evt, det):
        data = {}
        det.getData(evt)
        if det.evt.dat is None:
            return {}
        det.processFuncs()
        data = getUserData(det)
        data.update({"count": 1})  # to keep track of event count when combining in cube
        data = {det._name: data}
        return data

    def smd_default_detectors_engine(evt, det):
        return det.data(evt)

except ImportError:
    logger.error("Could not import smalldata_tools engine.")
