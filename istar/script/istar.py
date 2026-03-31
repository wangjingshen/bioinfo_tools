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

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, timer, run_with_single_thread

CONFIG = {
    "DEVICE": "cpu",  # or "cuda"
    "PixelSize": 0.5,
    "N_HVG": 1000,
    "EPOCHS": 400,
    "FilterSize": 8,
    "MinClusterSize": 20,
}


class IStar:
    def __init__(self, dir, image, spname, swap_pos, foreground_method, cluster_method, hard_radius, n_cluster, step):
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
        self.cluster_method = cluster_method
        self.hard_radius = hard_radius
        self.n_cluster = n_cluster
        self.spname = spname.rstrip("/") + "/"  #  spname+"/"  #istar prefix
        self.step = step.strip().split(',')


    @timer
    def input(self) -> None:
        '''
        generate_input
        '''
        mkdir(self.spname)
        execute_cmd((f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {ROOT}/gen_cnts_locs.R '
                f'--mtx {self.mtx} '
                f'--pos {self.pos} '
                f'--spname {self.spname} '
                f'--swap_pos {self.swap_pos}'
                ))
        if not os.path.exists("checkpoints"):
            execute_cmd(f'cp -r {ROOT}/../data/checkpoints/ .')

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
        execute_cmd(f'python {ROOT}/istar/rescale.py {self.spname}/ --image')
        execute_cmd(f'python {ROOT}/istar/preprocess.py {self.spname} --image')
        # extract histology features
        execute_cmd(f'python {ROOT}/istar/extract_features.py {self.spname} --device={CONFIG["DEVICE"]}')
        # select most highly variable genes to predict
        execute_cmd(f'python {ROOT}/istar/select_genes.py --n-top={CONFIG["N_HVG"]} {self.spname}/cnts.tsv {self.spname}/gene-names.txt')
        # rescale coordinates and spot radius
        execute_cmd(f'python {ROOT}/istar/rescale.py {self.spname} --locs --radius')

    
    @timer
    def get_mask(self) -> None:
        '''
        auto detect tissue mask
        '''
        execute_cmd(f'python {ROOT}/istar/get_mask_update.py {self.spname}/embeddings-hist.pickle {self.spname}/mask-small.png {self.foreground_method}')


    @timer
    def impute(self) -> None:
        '''
        predict super-resolution gene expression
        '''
        # train gene expression prediction model and predict at super-resolution
        if(self.hard_radius != None):
            # 
            with open(f"{self.spname}/radius.txt", "w") as f:
                f.write(self.hard_radius)

            execute_cmd(f'python {ROOT}/istar/impute.py {self.spname} --epochs={CONFIG["EPOCHS"]} --device={CONFIG["DEVICE"]}')
        else:
            # use 0 fix null 
            execute_cmd(f'python {ROOT}/istar/impute_update.py {self.spname} --epochs={CONFIG["EPOCHS"]} --device={CONFIG["DEVICE"]}')
        # visualize imputed gene expression
        execute_cmd(f'python {ROOT}/istar/plot_imputed.py {self.spname}')


    @timer
    def cluster(self) -> None:
        '''
        segment image by gene features
        '''
        #execute_cmd(f'python {ROOT}/istar/cluster.py --filter-size={CONFIG["FilterSize"]} --min-cluster-size={CONFIG["MinClusterSize"]} --n-clusters={CONFIG["N_CLUSTERS"]} --mask={self.spname}/mask-small.png {self.spname}/embeddings-gene.pickle {self.spname}/clusters-gene/')
        run_with_single_thread((f'python {ROOT}/istar/cluster.py ' 
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
        execute_cmd(f'python {ROOT}/istar/aggregate_imputed.py {self.spname}')
        execute_cmd(f'python {ROOT}/istar/reorganize_imputed.py {self.spname}')
        execute_cmd(f'python {ROOT}/istar/differential.py {self.spname}')
        # visualize spot-level gene expression data
        execute_cmd(f'python {ROOT}/istar/plot_spots.py {self.spname}')


    @timer
    def output(self) -> None:
        mkdir(f'{self.spname}/outs/')
        execute_cmd(f'cp -r  {self.spname}/clusters-gene/ {self.spname}/outs/')
        #execute_cmd(f'cd {self.spname}/outs/ && ln -s ../clusters-gene/ .')
        #execute_cmd(f'cd {self.spname}/outs/ && ln -s ../spots/ .')
        #execute_cmd(f'cd {self.spname}/outs/ && ln -s ../cnts-super-plots/ .')


    @timer
    def run(self) -> None:
        logger.info(f'{self.spname} start...')
        step_order = ['input', 'preprocess', 'get_mask', 'impute', 'cluster', 'downstream', 'output']
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
    parser.add_argument('--foreground_method', default = 'max', help='foreground method')
    parser.add_argument('--hard_radius', help='set hard radius')
    parser.add_argument('--cluster_method', default='km', help='cluster method')
    parser.add_argument('--n_cluster', default=10, help='number of cluster')
    parser.add_argument('--step', default='input,preprocess,get_mask,impute,cluster,downstream,output', help='comma-separated step')
    args = parser.parse_args()

    runner = IStar(args.dir, args.image, args.spname, args.swap_pos, args.foreground_method, 
                   args.cluster_method, args.hard_radius, args.n_cluster, args.step)
    runner.run()


if __name__ == '__main__':
    main()