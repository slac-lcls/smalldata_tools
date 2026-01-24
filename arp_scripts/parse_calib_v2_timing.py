#!/usr/bin/env python3
import argparse
import re
import statistics

def parse_time(msg):
    m = re.search(r'time=([0-9]+(?:\.[0-9]+)?)s', msg)
    return float(m.group(1)) if m else None

def add_value(stage_map, stage, val, rank):
    if val is None:
        return
    stage_map.setdefault(stage, []).append((val, rank))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "logfile",
        nargs="?",
        default="/sdf/data/lcls/ds/prj/prjlute22/scratch/for_chris_and_mona/smalldata_tools/arp_scripts/ds_count_calib_v2.log",
        help="Path to ds_count_calib_v2.log",
    )
    args = parser.parse_args()

    rank_re = re.compile(r'^\[DEBUG\]\s+Rank(\d+)\s+(.*)$')

    stage_values = {}
    with open(args.logfile, "r", errors="replace") as f:
        for line in f:
            m = rank_re.match(line)
            if not m:
                continue
            rank = int(m.group(1))
            msg = m.group(2).strip()
            if msg.startswith("Jungfrau calib cversion"):
                add_value(stage_values, "jungfrau_calib_cversion", parse_time(msg), rank)
                continue
            if msg.startswith("Jungfrau calib time"):
                add_value(stage_values, "jungfrau_calib_time", parse_time(msg), rank)
                continue

    print("units seconds")
    for key in ("jungfrau_calib_cversion", "jungfrau_calib_time"):
        entries = stage_values.get(key)
        if entries:
            vals = [v for v, _ in entries]
            avg = statistics.mean(vals)
            vmin = min(vals)
            vmax = max(vals)
            med = statistics.median(vals)
            std = statistics.pstdev(vals) if len(vals) > 1 else 0.0
            max_rank = next(r for v, r in entries if v == vmax)
            print(f"{key} {avg:.6f} {vmin:.6f} {vmax:.6f} {med:.6f} {std:.6f} {max_rank}")
        else:
            print(f"{key} NA NA NA NA NA NA")

if __name__ == "__main__":
    main()
