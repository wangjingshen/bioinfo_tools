import os
import pandas as pd
import time
import argparse
import subprocess
import time


class ScFEA():
    def __init__(self, args):
        '''
        version: 1.1.2
        '''
        self.rds = args.rds
        self.exp = args.exp
        self.model = args.model
        self.cell_type = args.cell_type
        self.group = args.group
        self.species = args.species
        self.n_threads = args.n_threads
        self.sc_imputation = args.sc_imputation
        self.moduleGene_file = args.moduleGene_file
        self.stoichiometry_matrix = args.stoichiometry_matrix
        self.cName_file = args.cName_file
        self.outdir = args.outdir
    
    def scFEA_input(self):
        cmd = f'Rscript /SGRNJ06/randd/USER/wangjingshen/project/scFEA1.1.2/script/scFEA_input.R \
                    --rds {self.rds} \
                    --exp {self.exp} \
                    --model {self.model} \
                    --cell_type {self.cell_type} \
                    --species {self.species} \
                    --group {self.group} \
                    --n_threads {self.n_threads} \
                    --sc_imputation {self.sc_imputation} \
                    --moduleGene_file {self.moduleGene_file} \
                    --stoichiometry_matrix {self.stoichiometry_matrix} \
                    --cName_file {self.cName_file} \
                    --outdir {self.outdir} '
        subprocess.call(cmd, shell=True)
    
    def scFEA_predication(self):
        #start = time.time() 
        df_mapfile = pd.read_csv(f"{self.outdir}/01.scFEA_input/mapfile.csv")     
        # run scFEA --
        for i in range(df_mapfile.shape[0]):
            print(df_mapfile.loc[i,'result_prefix'],"data load & process")
            cmd = f"python /SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/src/scFEA.py \
                        --data_dir /SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/data/ \
                        --test_file {df_mapfile.loc[i,'test_file']} \
                        --result_prefix {df_mapfile.loc[i,'result_prefix']} \
                        --moduleGene_file {df_mapfile.loc[i,'moduleGene_file']} \
                        --stoichiometry_matrix {df_mapfile.loc[i,'stoichiometry_matrix']} \
                        --cName_file {df_mapfile.loc[i,'cName_file']} \
                        --n_threads {df_mapfile.loc[i,'n_threads']} \
                        --sc_imputation {df_mapfile.loc[i,'sc_imputation']} \
                        --res_dir {self.outdir}/02.scFEA_predication/ "
            subprocess.call(cmd, shell=True)
        
        #end = time.time()
        #print("----\n") 
        #print("Running time:", end - start, "secs")

    def scFEA_visualization(self):
        cmd = f'Rscript /SGRNJ06/randd/USER/wangjingshen/project/scFEA1.1.2/script/scFEA_visualization.R \
                    --rds {self.rds} \
                    --model {self.model}  \
                    --cell_type {self.cell_type} \
                    --group {self.group} '
        subprocess.call(cmd, shell=True)

    def run(self):
        self.scFEA_input()
        self.scFEA_predication()
        #self.scFEA_visualization()


def main():
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--rds', type=str, help = 'seurat rds file', default = None)
    parser.add_argument('--exp', type=str, help = 'exp matrix file', default = None)
    parser.add_argument('--model', type=str, help = 'model, split or all', default = None)
    parser.add_argument('--cell_type', type=str, help = 'name of variable about cell type', default = None)
    parser.add_argument('--group', type=str, help = 'name of variable about group', default = None)
    parser.add_argument('--species', type=str, help = 'species', default = None)
    parser.add_argument('--n_threads', type=str, help = 'n_threads', default = 8)
    parser.add_argument('--sc_imputation', type=str, help = 'sc_imputation', default = "False")
    parser.add_argument('--moduleGene_file', type=str, help = 'moduleGene_file', default = None)
    parser.add_argument('--stoichiometry_matrix', type=str, help = 'stoichiometry_matrix', default = None)
    parser.add_argument('--cName_file', type=str, help = 'cName_file', default = None)
    parser.add_argument('--outdir', type=str, help = 'outdir', default = "./")
    args = parser.parse_args()

    runner = ScFEA(args)
    runner.run()


if __name__ == '__main__':
    main()