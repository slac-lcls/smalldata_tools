#!/usr/bin/env python
"""Print 1 if a producer config runs in integrating-detector mode for a run, else 0.

    python intg_probe.py <config> <run>

<config> is the producer config name as smd_producer.py takes it: the value of --config
(e.g. "qrix" for lcls2_producers/prod_config_qrix.py), or the module name itself. The
probe imports that module and prints 1 if get_intg(run)[0] is set, as smd_producer.py
does to decide on integrating mode. Any error prints 0 (and a short message on stderr),
so callers fall back to the non-integrating default.
"""

import contextlib
import importlib
import os
import sys


def main(argv):
    if len(argv) != 3:
        print(__doc__.strip(), file=sys.stderr)
        return "0"
    name, run = argv[1], argv[2]
    if not name.startswith("prod_config_"):
        name = f"prod_config_{name}"
    smd_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(smd_root, "lcls2_producers"))
    try:
        # Keep stdout for the answer only, whatever the config prints on import.
        with contextlib.redirect_stdout(sys.stderr):
            config = importlib.import_module(name)
            if not hasattr(config, "get_intg"):
                return "0"
            # smd_producer.py passes the run number as a string
            intg_main = config.get_intg(run)[0]
    except Exception as e:
        print(f"intg_probe: {name} run {run}: {type(e).__name__}: {e}", file=sys.stderr)
        return "0"
    return "1" if intg_main else "0"


if __name__ == "__main__":
    print(main(sys.argv))
