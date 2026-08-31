import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from contextlib import contextmanager
from larry import is_valid, in_filtered_list, larry


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def mkdir(dir) -> None:
    try:
        os.makedirs(dir, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")

def execute_cmd(command) -> None:
    logger.info(f"Executing: {command}")
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {command}, error: {e}")
        raise

@contextmanager
def tmp_chdir(path: Path):
    origin = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


class Larry:
    def __init__(self, fq_path, library_id, matrix, spname, step):        
        self.fq_path = fq_path
        self.library_id = library_id
        self.matrix = matrix
        self.spname = spname
        self.outdir = f'{self.spname}/outdir'
        self.step = step


    def input(self) -> None:
        '''
        generate_input
        '''
        mkdir(f'{self.spname}/input/matrix_10X')
        execute_cmd(f'cp {self.fq_path}/{self.library_id}*R1*f*gz {self.spname}/input/{self.library_id}_R1.fastq.gz')
        execute_cmd(f'cp {self.fq_path}/{self.library_id}*R2*f*gz {self.spname}/input/{self.library_id}_R2.fastq.gz')
        execute_cmd(f'cp {self.matrix} {self.spname}/input/matrix_10X/matrix_10X.tar')
        with tmp_chdir(f'{self.spname}/input/matrix_10X'):
            execute_cmd(f'tar -xvf matrix_10X.tar .')


    def larry_run(self) -> None:
        '''
        run larry
        '''
        mkdir(f'{self.outdir}')
        larry(self.spname, self.library_id, self.outdir)
    

    def run(self) -> None:
        logger.info(f'{self.spname} start...')
        step_order = ['input', 'larry_run']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()
        logger.info(f'{self.spname} completed.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fq_path', help='fq path', required=True)
    parser.add_argument('--library_id', help='library id', required=True)
    parser.add_argument('--matrix', help='matrix', required=True)
    parser.add_argument('--spname', help='spname', required=True)
    parser.add_argument('--step', default='input,larry_run', help='comma-separated step')
    args = parser.parse_args()

    runner = Larry(args.fq_path, args.library_id, args.matrix, args.spname, args.step)
    runner.run()


if __name__ == '__main__':
    main()