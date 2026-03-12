import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
from contextlib import contextmanager
import time

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
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