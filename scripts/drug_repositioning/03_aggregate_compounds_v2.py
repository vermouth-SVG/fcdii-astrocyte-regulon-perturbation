#!/usr/bin/env python3
from __future__ import annotations

import argparse

from drugrep_common_v2 import aggregate_compounds, diagnose_baseline, params_for


def main() -> None:
    parser = argparse.ArgumentParser(description="drug repurposing v2 aggregation")
    parser.add_argument("--run", choices=["baseline", "final"], default="baseline")
    args = parser.parse_args()
    params = params_for(args.run)
    aggregate_compounds(args.run, params)
    if args.run == "baseline":
        diagnose_baseline()


if __name__ == "__main__":
    main()
