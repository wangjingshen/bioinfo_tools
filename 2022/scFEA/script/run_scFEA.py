import os
import pandas as pd
import time
import argparse

all_start = time.time()

# args ----
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--rds', type=str, help = 'seurat rds file', default = None)
parser.add_argument('--model', type=str, help = 'model, split or all', default = None)
parser.add_argument('--cell_type', type=str, help = 'name of variable about cell type', default = None)
parser.add_argument('--group', type=str, help = 'name of variable about group', default = None)
parser.add_argument('--species', type=str, help = 'species', default = None)
parser.add_argument('--n_threads', type=str, help = 'n_threads', default = 8)
parser.add_argument('--sc_imputation', type=str, help = 'sc_imputation', default = "False")
parser.add_argument('--moduleGene_file', type=str, help = 'moduleGene_file', default = None)
parser.add_argument('--stoichiometry_matrix', type=str, help = 'stoichiometry_matrix', default = None)
parser.add_argument('--cName_file', type=str, help = 'cName_file', default = None)
args = parser.parse_args()

# 1 scFEA input ----
os.system("Rscript scFEA_input.R \
    --rds %s --model %s --cell_type %s --species %s --group %s --n_threads %s \
    --sc_imputation %s --moduleGene_file %s --stoichiometry_matrix %s --cName_file %s " \
    % (args.rds, args.model, args.cell_type, args.species, args.group, args.n_threads, 
       args.sc_imputation, args.moduleGene_file, args.stoichiometry_matrix, args.cName_file))


# 2 scFEA predication ----
print("\n*** 02.scFEA predication start running...")
start = time.time() 

df_mapfile = pd.read_csv("01.scFEA_input/mapfile.csv")     
# run scFEA --
for i in range(df_mapfile.shape[0]):
    print(df_mapfile.loc[i,'result_prefix'],"data load & process")
    os.system("python /SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/src/scFEA.py \
    --data_dir /SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/data/ \
    --test_file %s --result_prefix %s --moduleGene_file %s --stoichiometry_matrix %s \
    --cName_file %s --res_dir 02.scFEA_predication/ --n_threads %s --sc_imputation %s " \
    % (df_mapfile.loc[i,'test_file'],df_mapfile.loc[i,'result_prefix'],
    df_mapfile.loc[i,'moduleGene_file'], df_mapfile.loc[i,'stoichiometry_matrix'],
    df_mapfile.loc[i,'cName_file'], df_mapfile.loc[i,'n_threads'],
    df_mapfile.loc[i,'sc_imputation']))

end = time.time()
print("----\n") 
print("Running time:", end - start, "secs")   


# 3 scFEA visualization ----
os.system("Rscript scFEA_visualization.R \
    --rds %s --model %s --cell_type %s --group %s  " \
    % (args.rds, args.model, args.cell_type, args.group))

all_end = time.time()  
print("\n--------")
print("Total running time:", all_end - all_start, "secs. \n")