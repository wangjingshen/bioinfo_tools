#!/usr/bin/env python

import argparse
import subprocess


def update(q):
    '''
    cmd1: randd.q
    cmd2: ulimit -n 10000
    '''
    cmd1 = f"sed -i 's/sched_options -w n -cwd -V -l vf=100g,p=16/sched_options -q {q} -w n -cwd -V -l vf=200g,p=16/' sjm/sjm.job"
    cmd2 = f"sed -i 's/cmd source activate celescope2.2.0; celescope pathseq pathseq/cmd source activate celescope2.2.0; ulimit -n 10000; ulimit -n; celescope pathseq pathseq /' sjm/sjm.job"
    cmd3 = f"sed -i 's#cmd source activate celescope2.2.0#cmd source activate /SGRNJ/Public/Software/conda_env/celescope2.2.0#' sjm/sjm.job"
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    subprocess.check_call(cmd3, shell=True)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--q', required=True)

    args = parser.parse_args()
    update(args.q)
    print(f'queue: {args.q}')
    print("sjm update done.")