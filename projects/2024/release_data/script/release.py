import pandas as pd
import argparse
import subprocess

class Release:
    def __init__(self, args):
        self.mapfile = args.mapfile
        self.analysis_path = args.analysis_path

    def cp(self):
        '''
        L240925063      /SGRNJ06/DATA04/limsdownload/24_10/2024_10_12/P24041601/N-11577100/2024-10-12-61        N_11577100_RNA
        '''
        mapfile = pd.read_csv(self.mapfile, header=None, sep="\t")
        
        for i in range(mapfile.shape[0]):
            cmd1 = f'mkdir -p release/fastq/{mapfile.iloc[i,2]}'
            cmd2 = f'cp -r {mapfile.iloc[i,1]}/{mapfile.iloc[i,0]}*gz  release/fastq/{mapfile.iloc[i,2]}'
            subprocess.check_call(cmd1, shell = True)
            subprocess.check_call(cmd2, shell = True)
            
            if(self.analysis_path != None):
                cmd3 = f'mkdir -p release/mtx/{mapfile.iloc[i,2]}'
                cmd4 = f'cp -r {self.analysis_path}/{mapfile.iloc[i,2]}/outs/filtered  release/mtx/{mapfile.iloc[i,2]}'
                subprocess.check_call(cmd3, shell = True)
                subprocess.check_call(cmd4, shell = True)


    def run(self):
        self.cp()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile')
    parsers.add_argument('--analysis_path', help='analysis_path')

    args = parsers.parse_args()
    runner = Release(args) 
    runner.run()

if __name__ == '__main__':
    main()
