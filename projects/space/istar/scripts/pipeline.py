import os
os.environ['OMP_NUM_THREADS'] = '32'
os.environ['MKL_NUM_THREADS'] = '32'
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
#os.environ['NUMEXPR_NUM_THREADS'] = '1'
import sys
import pandas as pd
import argparse
import glob
import subprocess
import logging
import json
import psutil

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def add_root():
    root = Path(__file__).resolve().parent
    while not (root / "utils").exists():
        parent = root.parent
        if parent == root:
            raise FileNotFoundError("utils not found！")
        root = parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

bioinfo_root = add_root()
pipeline_root = Path(__file__).resolve().parents[1]

from utils.utils import mkdir, logger, execute_cmd, timer, run_with_single_thread
from istarSpots import istar2spots

CONFIG = {
    "DEVICE": "cpu",  # or "cuda"
    "PixelSize": 0.5,
    "N_HVG": 1000,
    "EPOCHS": 400,
    "FilterSize": 8,
    "MinClusterSize": 20,
}


class IStar:
    def __init__(self, dir, image, spname, swap_pos, foreground_method, foreground_cluster_method, cluster_method, hard_radius, n_cluster, k, distance_thresh, clip, step):
        self.dir = dir
        filtered_dir = os.path.join(dir, "outs", "filtered")
        pos = os.path.join(dir, "outs", "spatial", "positions_list.csv")

        if not os.path.exists(filtered_dir):
            raise FileNotFoundError(f"[IStar] Required directory not found: {filtered_dir}")
        if not os.path.isfile(image):
            raise FileNotFoundError(f"[IStar] Image file not found: {image}")
        if not os.path.isfile(pos):
            raise FileNotFoundError(f"[IStar] Position file not found: {pos}")
        
        self.mtx = filtered_dir
        self.image = image
        self.pos = pos

        self.swap_pos = swap_pos
        self.foreground_method = foreground_method
        self.foreground_cluster_method = foreground_cluster_method
        self.cluster_method = cluster_method
        self.hard_radius = hard_radius
        self.n_cluster = n_cluster
        self.spname = spname.rstrip("/") + "/"  #  spname+"/"  #istar prefix
        self.k = k
        self.distance_thresh = distance_thresh
        self.clip = clip
        self.step = step.strip().split(',')


    @timer
    def input(self) -> None:
        '''
        generate_input
        '''
        mkdir(self.spname)
        execute_cmd((f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {pipeline_root}/scripts//gen_cnts_locs.R '
                f'--mtx {self.mtx} '
                f'--pos {self.pos} '
                f'--spname {self.spname} '
                f'--swap_pos {self.swap_pos}'
                ))
        if not os.path.exists("checkpoints"):
            execute_cmd(f'cp -r {pipeline_root}/scripts//../data/checkpoints/ .')

        with open(f'{self.dir}/outs/spatial/scalefactors_json.json') as f:
            sf = json.load(f)
        PixelSizeRaw = 8000/2000 * float(sf['tissue_hires_scalef'])  #https://github.com/daviddaiweizhang/istar
        RADIUS_RAW = 0.5 * float(sf['spot_diameter_fullres'])

        with open(f"{self.spname}/radius-raw.txt", "w") as f:
            f.write(str(RADIUS_RAW))
        with open(f"{self.spname}/pixel-size-raw.txt", "w") as f:
            f.write(str(PixelSizeRaw))
        with open(f"{self.spname}/pixel-size.txt", "w") as f:
            f.write(str(CONFIG["PixelSize"]))

        execute_cmd(f'cp {self.image} {self.spname}/he-raw.png')


    @timer
    def preprocess(self) -> None:
        '''
        preprocess_image
        '''
        execute_cmd(f'python {pipeline_root}/scripts//istar/rescale.py {self.spname}/ --image')
        execute_cmd(f'python {pipeline_root}/scripts//istar/preprocess.py {self.spname} --image')
        # extract histology features
        execute_cmd(f'python {pipeline_root}/scripts//istar/extract_features.py {self.spname} --device={CONFIG["DEVICE"]}')
        # select most highly variable genes to predict
        execute_cmd(f'python {pipeline_root}/scripts//istar/select_genes.py --n-top={CONFIG["N_HVG"]} {self.spname}/cnts.tsv {self.spname}/gene-names.txt')
        # rescale coordinates and spot radius
        execute_cmd(f'python {pipeline_root}/scripts//istar/rescale.py {self.spname} --locs --radius')

    
    @timer
    def getMask(self) -> None:
        '''
        istar: auto detect tissue mask
        in_tissue: get mask from positions_list
        '''
        if(self.foreground_method == "istar"):
            execute_cmd(f'python {pipeline_root}/scripts//istar/get_mask_update.py {self.spname}/embeddings-hist.pickle {self.spname}/mask-small.png {self.foreground_cluster_method}')
        
        if(self.foreground_method == "in_tissue"):
            execute_cmd(f'python {pipeline_root}/scripts//istar/get_mask_in_tissue.py {self.pos} {self.spname}/mask-small.png')


    @timer
    def impute(self) -> None:
        '''
        predict super-resolution gene expression
        '''
        if(self.hard_radius != None):
            with open(f"{self.spname}/radius.txt", "w") as f:
                f.write(self.hard_radius)

        # train gene expression prediction model and predict at super-resolution
        execute_cmd(f'python {pipeline_root}/scripts//istar/impute_update.py {self.spname} --epochs={CONFIG["EPOCHS"]} --device={CONFIG["DEVICE"]}')
        # visualize imputed gene expression
        execute_cmd(f'python {pipeline_root}/scripts//istar/plot_imputed.py {self.spname}')


    @timer
    def cluster(self) -> None:
        '''
        segment image by gene features
        '''
        #execute_cmd(f'python {pipeline_root}/scripts//istar/cluster.py --filter-size={CONFIG["FilterSize"]} --min-cluster-size={CONFIG["MinClusterSize"]} --n-clusters={CONFIG["N_CLUSTERS"]} --mask={self.spname}/mask-small.png {self.spname}/embeddings-gene.pickle {self.spname}/clusters-gene/')
        run_with_single_thread((f'python {pipeline_root}/scripts//istar/cluster.py ' 
                          f'--filter-size={CONFIG["FilterSize"]} '
                          f'--min-cluster-size={CONFIG["MinClusterSize"]} '
                          f'--n-clusters={self.n_cluster} '
                          f'--mask={self.spname}/mask-small.png '
                          f'--method={self.cluster_method} '
                          f'{self.spname}/embeddings-gene.pickle '
                          f'{self.spname}/clusters-gene/'))
    

    @timer
    def downstream(self) -> None:
        '''
        differential & visualize
        '''
        # differential analysis by clusters
        execute_cmd(f'python {pipeline_root}/scripts//istar/aggregate_imputed.py {self.spname}')
        execute_cmd(f'python {pipeline_root}/scripts//istar/reorganize_imputed.py {self.spname}')
        execute_cmd(f'python {pipeline_root}/scripts//istar/differential.py {self.spname}')
        # visualize spot-level gene expression data
        execute_cmd(f'python {pipeline_root}/scripts//istar/plot_spots.py {self.spname}')


    @timer
    def output(self) -> None:
        mkdir(f'{self.spname}/outs/')
        execute_cmd(f'cp -r  {self.spname}/clusters-gene/ {self.spname}/outs/istar')


    @timer
    def istarSpots(self) -> None:
        '''
        {self.outdir}/istar2spots.png
        {self.outdir}/istar2spots_df.csv
        '''
        logger.info(f"[{self.spname}] Running ln_outs step.")

        self.space_input = f'{self.spname}/space_input/'
        mkdir(self.space_input)
        execute_cmd((f'cp {self.dir}/outs/filtered_feature_bc_matrix.h5 {self.space_input}'))
        execute_cmd((f'cp -r {self.dir}/outs/spatial {self.space_input}'))
        execute_cmd((f'mv {self.space_input}/spatial/positions_list.csv {self.space_input}/spatial/tissue_positions_list.csv'))
        self.space_pos = f'{self.space_input}/spatial/tissue_positions_list.csv'

        self.istar_labels = f'{self.spname}/outs/istar/labels.pickle'
        self.istarSpotsOuts = f'{self.spname}/outs/istarSpots'
        mkdir(self.istarSpotsOuts)
        istar2spots(self.istar_labels, self.space_pos, self.istarSpotsOuts, self.foreground_method, self.k, self.distance_thresh, self.clip)

        execute_cmd(f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {pipeline_root}/scripts/seurat_plot.R '
                    f'--space_input {self.space_input} '
                    f'--istar2spots {self.istarSpotsOuts}/istar2spots_df.csv '
                    f'--outdir {self.istarSpotsOuts} ')

    @timer
    def orgdir(self) -> None:
        mkdir(f'run/')
        execute_cmd(f'mv {self.spname}/* run ')
        execute_cmd(f'mv run {self.spname} ')
        execute_cmd(f'mv {self.spname}/run/outs  {self.spname}/')
        execute_cmd(f'rm -rf {self.spname}/run/space_input ')
        execute_cmd(f'rm -rf checkpoints ')


    @timer
    def run(self) -> None:
        logger.info(f'{self.spname} start...')
        step_order = ['input', 'preprocess', 'getMask', 'impute', 'cluster', 'downstream', 'output', 'istarSpots', 'orgdir']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()
        logger.info(f'{self.spname} completed.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', help='celescope dir', required=True)
    parser.add_argument('--image', help='image', required=True)
    parser.add_argument('--spname', help='spname', required=True)
    parser.add_argument('--swap_pos', default = "T", help='swap postion')
    parser.add_argument('--foreground_method', default = 'in_tissue', help='foreground method, in_tissue or istar')
    parser.add_argument('--foreground_cluster_method', default = 'max', help='foreground cluster method in istar, max or min')
    parser.add_argument('--hard_radius', help='set hard radius')
    parser.add_argument('--cluster_method', default='km', help='cluster method')
    parser.add_argument('--n_cluster', default=10, help='number of cluster')
    parser.add_argument('--k', help='k istar pixels', default = 3, type = int)
    parser.add_argument('--distance_thresh', help='thresh of distance between spots and istar pixels', default=200, type = int)
    parser.add_argument('--clip', help='clip spots from istar', action='store_true') 
    parser.add_argument('--step', default='input,preprocess,getMask,impute,cluster,downstream,output,istarSpots,orgdir', help='comma-separated step')
    args = parser.parse_args()

    runner = IStar(args.dir, args.image, args.spname, args.swap_pos, args.foreground_method, args.foreground_cluster_method,
                   args.cluster_method, args.hard_radius, args.n_cluster, args.k, args.distance_thresh, args.clip, args.step)
    runner.run()


if __name__ == '__main__':
    main()