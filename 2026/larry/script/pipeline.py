import os
import sys
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))

from larry import is_valid, in_filtered_list, larry
from utils.utils import mkdir, logger, execute_cmd, timer, tmp_chdir


class Larry:
    def __init__(self, fq_path, library_id, matrix, spname, step):        
        self.fq_path = fq_path
        self.library_id = library_id
        self.matrix = matrix
        self.spname = spname
        self.outdir = f'{self.spname}/outdir'
        self.step = step


    @timer
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


    @timer
    def larry_run(self) -> None:
        '''
        run larry
        '''
        mkdir(f'{self.outdir}')
        larry(self.spname, self.library_id, self.outdir)
    

    @timer
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