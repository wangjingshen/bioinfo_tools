import os
import argparse
import subprocess
import pandas as pd


class Merge_fastq:
    def __init__(self, args):
        self.mapfile = args.mapfile
        self.pipeline = args.pipeline

        if not os.path.exists("cmd"):
            os.system(f"mkdir -p cmd")
        if not os.path.exists("raw_data"):
            os.system(f"mkdir -p raw_data")

    def generate_merge_cmd(self):
        df_mapfile = pd.read_csv(self.mapfile, header= 0, sep="\t")        
        def generate_cmd(reads):
            cmd = df_mapfile.iloc[0,1] + "/" + df_mapfile.iloc[0,0] + "_" + reads + ".fastq.gz "
            for i in range(len(df_mapfile)-1) :
                cmd = cmd + df_mapfile.iloc[i+1,1] + "/" + df_mapfile.iloc[i+1,0] + "_" + reads + ".fastq.gz "
            cmd = "cat "+ cmd + "> raw_data/" + df_mapfile.iloc[0,0] + "_" + reads + ".fastq.gz"
            return(cmd)

        cmd1 = generate_cmd("R1")
        cmd2 = generate_cmd("R2")
        with open(f"cmd/merge_R1.sh", 'w') as f:
            f.write(cmd1)
        with open(f"cmd/merge_R2.sh", 'w') as f:
            f.write(cmd2)
    
    def run_merge_cmd(self):
        cmd1 = f"bash cmd/merge_R1.sh"
        cmd2 = f"bash cmd/merge_R2.sh"
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
    
    def generate_mapfile(self):
        df_mapfile = pd.read_csv(self.mapfile, header= 0, sep="\t").head(1)   # first row 
        path = os.getcwd()    
        df_mapfile.iloc[0,1] = f"{path}/raw_data/"     # absolute path
        df_mapfile.to_csv("celescope_mapfile", sep="\t", header=None, index=None)
    
    def generate_celescope_cmd(self):
        if(self.pipeline == "HBV_fj"):
            cmd1 = f"source activate celescope1.5.1b0"
            cmd2 = f"multi_capture_virus\
                    --mapfile celescope_mapfile \
                    --virus_genome /SGRNJ03/randd/test_rd/dxh/data/HBV_genome/ \
                    --thread 4\
                    --not_consensus \
                    --outFilterMatchNmin 80 \
                    --umi_threshold otsu \
                    --allowNoPolyT \
                    --gtf /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-10-18_gtf/HBV.gtf"
        
        with open(f"run_celescope.sh", 'w') as f:
            f.write(cmd1 + "\n")
            f.write(cmd2)
    
    def run_celescope_cmd(self):
        cmd = f"bash run_celescope.sh"
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.generate_merge_cmd()
        self.run_merge_cmd()
        self.generate_mapfile()
        self.generate_celescope_cmd()
        self.run_celescope_cmd()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile', required=True)
    parsers.add_argument('--pipeline', help='pipeline', required=True)


    args = parsers.parse_args()
    runner = Merge_fastq(args) 
    runner.run()

if __name__ == '__main__':
    main()





