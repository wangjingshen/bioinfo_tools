import argparse
import os
import scanpy as sc
import pandas as pd

def sum_pos(x):
    return(sum(x>0))

def get_len(x):
    return(len(x.to_list()))

def HBV_stat(h5ad, virus_UMI, outdir):
    if(os.path.exists(f'{outdir}')):
        os.system(f'mkdir -p {outdir}')
    sc.settings.figdir = outdir

    h5ad = sc.read_h5ad(h5ad)
    virus_UMI = pd.read_csv(virus_UMI,index_col=0).fillna(0)

    # check colnames and index
    if (sum(h5ad.obs_names == virus_UMI.index) == h5ad.n_obs):
        h5ad.obs['HBV_UMI'] = virus_UMI.sum_UMI
        h5ad.obs['HBV_status'] = h5ad.obs['HBV_UMI'].apply(lambda x: 'HBV+' if x>0 else 'HBV-').astype('category')

        sc.pl.umap(h5ad, color='HBV_status', palette=["red","lightgrey"], show = False, save="_HBV_status")

        stat_df = h5ad.obs.groupby('cluster').agg({
            'HBV_UMI': sum_pos,
            'HBV_status': get_len
        }).rename(columns={"HBV_UMI": "HBV_postive_cell_number",
                           "HBV_status": "cell_number"}).reset_index()
        stat_df['HBV_positive_percent'] = round(stat_df.HBV_postive_cell_number/stat_df.cell_number*100, 2)

        stat_df.to_csv(f'{outdir}/HBV_stat.tsv', sep="\t", index = False)
    else:
        print("zl and fj no match.")

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--h5ad', help='h5ad', required = True)
    parsers.add_argument('--virus_UMI', help="virus_UMI", required = True)
    parsers.add_argument('--outdir', help="outdir", required = True)

    args = parsers.parse_args()
    HBV_stat(args.h5ad, args.virus_UMI, args.outdir)

if __name__ == '__main__':
    main()