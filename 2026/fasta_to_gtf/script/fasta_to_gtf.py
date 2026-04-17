#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
import time

from bam_to_gtf import bam_to_bed,bed_to_gtf

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, tmp_chdir


class FastaToGtf:
    def __init__(self, genome:Path, fasta:Path, prefix:str, star_path:Path, threads:int=4, force:bool=False):
        self.genome = genome.resolve()
        self.fasta = fasta.resolve()
        self.prefix = prefix
        self.star_path = star_path
        self.threads = threads
        self.force = force
        self.star_out = Path("star")
        self.gtf_out = Path(f'{prefix}.gtf')


    def star_align(self) -> None:
        """STAR to get BAM"""
        bam_file = self.star_out / f"{self.prefix}Aligned.sortedByCoord.out.bam"
        if bam_file.exists() and not self.force:
            logger.info(f"BAM already exists: {bam_file}  (--force to overwrite)")
            return

        mkdir(self.star_out)
        cmd = (
            f"{self.star_path} "
            f"--runThreadN {self.threads} "
            f"--genomeDir {self.genome} "
            f"--readFilesIn {self.fasta} "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outFileNamePrefix {self.prefix}"
        )
        logger.info("Running STAR...")
        with tmp_chdir(self.star_out):
            execute_cmd(cmd)
        self.bam = f'{self.star_out}/{self.prefix}Aligned.sortedByCoord.out.bam'


    def bam_to_gtf(self) -> None:
        """BAM -> bam -> gtf"""
        bed = Path(f"{self.prefix}.bed")
        bam_to_bed(self.bam, bed)

        gtf = Path(f"{self.prefix}.gtf")
        bed_to_gtf(bed, gtf, self.prefix)


    def update(self, bed: Path) -> None:
        """update bed -> GTF"""
        gtf = Path(f"{self.prefix}.gtf")
        bed_to_gtf(bed, gtf)


    def run(self) -> None:
        t0 = time.time()
        self.star_align()
        self.bam_to_gtf()
        logger.info(f"Done. ({time.time() - t0:.1f}s)")
    

def update_bed_to_gtf(prefix, bed: Path) -> None:
    """update bed -> GTF"""
    t0 = time.time()
    gtf = Path(f"{prefix}_update.gtf")
    bed_to_gtf(bed, gtf, prefix)
    logger.info(f"Done. ({time.time() - t0:.1f}s)")


def main():
    parser = argparse.ArgumentParser(description="FASTA → STAR → BAM → BED → GTF")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common_args(p):
        p.add_argument("--prefix", required=True, help="prefix")

    p_run = sub.add_parser("run", help="run full pipeline")
    add_common_args(p_run)
    p_run.add_argument("--genome", required=True, type=Path, help="genome path")
    p_run.add_argument("--fasta",  required=True, type=Path, help="fasta")
    p_run.add_argument("--star",   default="/Public/Software/miniconda2/envs/old/bin/STAR",
        help="STAR path, /SGRNJ/Public/Software/conda_env/celescope2.1.0/bin/STAR-avx2,\
                         /Public/Software/miniconda2/envs/old/bin/STAR")
    p_run.add_argument("--threads", default=4, type=int, help="threads")
    p_run.add_argument("--force",  action="store_true", help="force get bam")

    p_update = sub.add_parser("update", help="use updated bed to GTF")
    add_common_args(p_update)
    p_update.add_argument("--bed", required=True, type=Path, help="updated bed path")

    args = parser.parse_args()

        
    if args.command == "run":
        star_path = Path(args.star).expanduser()
        if not star_path.is_file():
            raise FileNotFoundError(f"STAR not found: {star_path}")
        runner = FastaToGtf(args.genome, args.fasta, args.prefix, star_path, args.threads, args.force)
        runner.run()
    else:  # if args.command == "update":
        update_bed_to_gtf(args.prefix, args.bed)


if __name__ == "__main__":
    main()