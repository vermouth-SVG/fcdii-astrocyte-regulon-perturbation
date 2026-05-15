#!/usr/bin/env python3
from __future__ import annotations

import argparse

from drugrep_common_v2 import params_for, run_enrichment


def main() -> None:
    parser = argparse.ArgumentParser(description="drug repurposing v2 DSigDB/Enrichr")
    parser.add_argument("--run", choices=["baseline", "final"], default="baseline")
    args = parser.parse_args()
    params = params_for(args.run)
    run_enrichment(args.run, params)


if __name__ == "__main__":
    main()
