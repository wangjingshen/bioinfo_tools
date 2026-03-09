import os
os.environ['OMP_NUM_THREADS'] = '32' 
os.environ['MKL_NUM_THREADS'] = '32'

import sys
import argparse
import scanpy as sc

from pathlib import Path

from h5toh5ad import h5toh5ad
from cut_visium_spots import cut_visium_spots
from encode import get_top_k_genes, encode_mtx, encode_image, finetune

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, timer

class Loki_preprocess:
    def __init__(self, dir: Path, spname: str, hk_genes: str, sc_h5ad: str, step: str):
        self.dir = dir
        self.spname = spname
        self.space_input = f'{self.spname}/space_input/'
        self.hk_genes = hk_genes
        self.sc_h5ad = sc_h5ad
        self.loki_input = f'{self.spname}/loki_input/'
        self.loki_input_finetune = f'{self.spname}/loki_input_finetune/'
        self.step = step.strip().split(',')

    def input(self) -> None:
        '''
        generate_input
        '''
        logger.info(f"[{self.spname}] Running space_input step.")
        mkdir(self.space_input)
        mkdir(self.loki_input)

        execute_cmd((f'cp {self.dir}/outs/filtered_feature_bc_matrix.h5 {self.space_input}'))
        execute_cmd((f'cp -r {self.dir}/outs/spatial {self.space_input}'))
        execute_cmd((f'mv {self.space_input}/spatial/positions_list.csv {self.space_input}/spatial/tissue_positions_list.csv'))

        invalid_spots = cut_visium_spots(self.space_input, f'{self.loki_input}/spots_images/')
        logger.info(f"[{self.spname}] invalid spots:{invalid_spots}.")

        logger.info(f"[{self.spname}] Running h5toh5ad step.")
        self.h5ad = h5toh5ad(self.space_input, invalid_spots)
        logger.info(f"[{self.spname}] space_input done.")
            
    def loki_encode(self) -> None:
        '''
        loki_encode
        '''
        logger.info(f"[{self.spname}] Running loki_encode step.")
        self.top_k_genes = get_top_k_genes(self.h5ad, self.hk_genes)
        encode_mtx(self.top_k_genes, self.loki_input, "space")
        encode_image(self.space_input, f'{self.space_input}/image_coord.csv', self.loki_input)

        if self.sc_h5ad is not None:
            sc_h5ad = sc.read_h5ad(self.sc_h5ad)
            sc_h5ad.var_names_make_unique()
            top_k_genes_sc = get_top_k_genes(sc_h5ad, self.hk_genes)
            encode_mtx(top_k_genes_sc, self.loki_input, "sc")

        logger.info(f"[{self.spname}] loki_encode done.")
    
    def loki_finetune_encode(self) -> None:
        '''
        loki_finetune_encode
        '''
        logger.info(f"[{self.spname}] Running loki_finetune_encode step.")
        mkdir(self.loki_input_finetune)
        finetune(self.top_k_genes, f'{self.loki_input}/spots_images/', self.loki_input_finetune, self.spname)
        finetune_model = f'logs/finetune_{self.spname}/checkpoints/epoch_10.pt'
        invalid_spots = cut_visium_spots(self.space_input, f'{self.loki_input_finetune}/spots_images/')
        encode_mtx(self.top_k_genes, self.loki_input_finetune, "space", model_path = finetune_model)
        encode_image(self.space_input, f'{self.space_input}/image_coord.csv', self.loki_input_finetune, model_path = finetune_model)

        logger.info(f"[{self.spname}] loki_finetune_encode done.")

    @timer
    def run(self) -> None:
        step_order = ['input','loki_encode','loki_finetune_encode']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', help='celescope sapce dir', required=True)
    parser.add_argument('--spname', help='spname', required=True)
    parser.add_argument('--hk_genes', help='housekeeping_genes', required=True)
    parser.add_argument('--sc_h5ad', help='sc h5ad', default=None)
    parser.add_argument('--step', default='input,loki_encode,loki_finetune_encode', help='comma-separated step')
    args = parser.parse_args()

    runner = Loki_preprocess(args.dir, args.spname, args.hk_genes, args.sc_h5ad, args.step)
    runner.run()


if __name__ == '__main__':
    main()