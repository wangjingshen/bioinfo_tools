import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
from contextlib import contextmanager
import time

def add_root(levels_up=5):
    root = Path(__file__).resolve()
    for _ in range(levels_up):
        root = root.parent
    if not (root / "utils").exists():
        raise FileNotFoundError(f"utils not found in {root}.")
    sys.path.insert(0, str(root))
    return(root)

root_path = add_root(5)  # Top 5 parent directories of current script (bioinfo_tools)
script_path = Path(__file__).resolve().parent

from utils.utils import mkdir, logger, execute_cmd


def bam_to_bed(bam, bed) -> None:
    """BAM → bed12"""
    logger.info("BAM → BED...")
    cmd = f"bedtools bamtobed -i {bam} -bed12 > {bed}"
    execute_cmd(cmd)

def bed_to_gtf(bed: Path, gtf: Path, name) -> None:
    """bed → genePred → GTF"""
    logger.info("BED → GTF...")
    cmd = (
        f"bedToGenePred {bed} /dev/stdout | "
        f"genePredToGtf file /dev/stdin tmp && "
        f"sed -i 's|/dev/stdin|{name}_genome|g' tmp && "
        f"awk -F'\\t' -v OFS='\\t' '{{if ($3==\"transcript\") $3=\"gene\"; print}}' tmp > {gtf} && "
        f"rm -f tmp"
    )
    execute_cmd(cmd)