import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detectors = []


def get_intg(run):
    """
    Returns
    -------
    intg_main (str):  This detector ill be passed to the psana datasource. It should be the SLOWEST of
                all integrating detectors in the data
    intg_addl (list of str): The detectors in this list will be analyzed as integrating detectors. It is
                             important that the readout frequency of these detectors is commensurate and
                             in-phase with intg_main.
    """
    run = int(run)
    intg_main = None
    intg_addl = []
    if run > 0:
        intg_main = ""
        intg_addl = []
    return intg_main, intg_addl


slow_detectors = []  # NOT IMPLEMENTED


##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
# epicsPV = ['las_fs14_controller_time']
# epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
# epicsOncePV = [('GDET:FEE1:241:ENRC', "MyTest"), 'GDET:FEE1:242:ENRC', "FOO:BAR:BAZ"]
epicsPV = []
epicsOncePV = []


##########################################################
# psplot config
##########################################################

import psplot


def get_psplot_configs(run):
    configs = {}

    if run > 0:
        andor_vls = {
            "callback": psplot.SpectrumScan,
            "data_fields": {
                "data": "andor_vls/unaligned_hitfinder",
                "norm": "andor_vls/unaligned_norm_det_rix_fim1_sum_wfintegrate",
                "count": "andor_vls/unaligned_norm_count",
                "scan": "andor_vls/unaligned_norm_scan_sum_mono_ev",
            },
            "bins": np.linspace(520, 540, 30),
            "spectrum_range": (900, 1050),
            "lineouts_idx": (50, 70, 100),
        }
        configs["andor_vls"] = andor_vls
    return configs
