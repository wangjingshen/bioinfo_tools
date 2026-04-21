import os
import subprocess
import argparse

def md5sum(path):
    '''
    # v0 ----
    In release_data.py, copy commands are run in the background.
    Some data may be not copied completely, so a separate script is written for the md5sum command.
    run in the work dir of release_data
    python /SGRNJ06/randd/USER/wangjingshen/script/release_data/md5sum.py

    # v1 ----
    not need
    '''
    os.chdir(path)   # chdir for write md5sum.txt in data/
    cmd1 = f"find fastq/ -type l -print0 | xargs -0 md5sum > md5sum.txt"
    cmd2 = f"zip md5sum.txt.zip md5sum.txt"    # zip for oss
    subprocess.check_call(cmd1, shell = True)
    subprocess.check_call(cmd2, shell = True)

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--path', help='path', required=True)

    args = parsers.parse_args()
    md5sum(args.path)

if __name__ == '__main__':
    main()