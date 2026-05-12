import os
import sys
import argparse
from pathlib import Path

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