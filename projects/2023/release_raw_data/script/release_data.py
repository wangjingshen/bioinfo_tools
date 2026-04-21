import os
import subprocess
import argparse


class Release_data():
    def __init__(self, args):
        self.mapfile = args.mapfile
        self.project_id = args.project_id
        self.subset = args.subset
        self.mode = args.mode
        self.rm_version = args.rm_version
        self.outdir = args.outdir
        self.generate_mapfile = args.generate_mapfile
        self.generate_md5sum = args.generate_md5sum
        self.time = args.time

    def generate_cmd(self):
        cmd = f"Rscript /SGRNJ06/randd/USER/wangjingshen/script/release_raw_data/generate_mkdir_ln_cmd.R \
                    --mapfile {self.mapfile} \
                    --project_id {self.project_id} \
                    --subset {self.subset} \
                    --mode {self.mode} \
                    --rm_version {self.rm_version} \
                    --outdir {self.outdir} \
                    --generate_mapfile {self.generate_mapfile} \
                    --time {self.time}"
        subprocess.check_call(cmd, shell = True)
    
    def run_mkdir(self):
        if(self.mode == "zl"):
            cmd = f"bash cmd/zl_mkdir.sh"
        if(self.mode == "fj"):
            cmd = f"bash cmd/fj_mkdir.sh"
        subprocess.check_call(cmd, shell = True)
    
    def run_ln(self):
        if(self.mode == "zl"):
            cmd = f"bash cmd/zl_ln.sh"
        if(self.mode == "fj"):
            cmd = f"bash cmd/fj_ln.sh"
        subprocess.check_call(cmd, shell = True)
    
    def run_md5sum(self):
        cmd = f'python /SGRNJ06/randd/USER/wangjingshen/script/release_raw_data/md5sum.py --path {self.outdir}'
        subprocess.check_call(cmd, shell = True)
    
    def check_md5sum(self):
        cmd1 = f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/release_raw_data/check_md5sum.R --path {self.outdir}'
        cmd2 = f'rm {self.outdir}/md5sum.txt'
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)

    def run(self):
        self.generate_cmd()
        self.run_mkdir()
        self.run_ln()

        if(self.generate_md5sum != 'no'):
            self.run_md5sum()
            self.check_md5sum()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile, split by ,', required=True)
    parsers.add_argument('--project_id', help='project id')
    parsers.add_argument('--subset', help='subset, split by , ')
    parsers.add_argument('--mode', help='mode, split by , ')
    parsers.add_argument('--rm_version', help='T or F')
    parsers.add_argument('--outdir', help='outdir')
    parsers.add_argument('--generate_mapfile', help='generate_mapfile',default='no')
    parsers.add_argument('--generate_md5sum', help='generate_md5sum',default='no')
    parsers.add_argument('--time', help='time')

    args = parsers.parse_args()
    runner = Release_data(args) 
    runner.run()

if __name__ == '__main__':
    main()
        