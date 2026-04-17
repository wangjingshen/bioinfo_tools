import os
import pandas as pd
import argparse
import json
import glob


class Json_to_df():
    def __init__(self, args):
        self.prefix = args.prefix
        if(args.mod == "genome"):
            self.json_file = glob.glob(f'{args.cele_path}/*06.filter_virus*/*filtered_read_count.json')[0]
        if(args.mod == "split"):
            self.json_file = glob.glob(f'{args.cele_path}/*06.filter_virus*/*corrected_read_count.json')[0]
        self.outdir = args.outdir
        if not os.path.exists(self.outdir}:
            os.system(f'mkdir -p {self.outdir}')
        
    
    def run_json_to_df(self):
        f = open(f'{self.json_file}','r')
        json_file = json.loads(f.read())
        f.close()
        df_file = pd.json_normalize(json_file)
        df_file.T.to_csv(f'{self.outdir}/{self.prefix}_.tsv"', sep="\t", header=None)

    def run(self):
        self.run_json_to_df()
        

def main():
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--cele_path', help = 'celescope path', required = True)
    parser.add_argument('--prefix', help = 'prefix', required = True)
    parser.add_argument('--mod', help = 'mod', required = True)
    parser.add_argument('--outdir', help = 'outdir', required = True)
    args = parser.parse_args()
    runner = Json_to_df(args)
    runner.run()


if __name__ == '__main__':
    main()