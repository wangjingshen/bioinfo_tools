import os
import sys
import torch
import subprocess
import open_clip
import pandas as pd
import numpy as np
import scanpy as sc

import loki.utils
import loki.align
import loki.preprocess
import loki.utils
import loki.plot

import matplotlib.pyplot as plt

from PIL import Image
from pathlib import Path
from matplotlib.colors import rgb2hex


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = f'{ROOT}/data/checkpoint.pt'
MODEL_NAME = 'coca_ViT-L-14'
DEVICE = 'cpu'


def get_top_k_genes(h5ad, housekeeping_genes):
    '''
    get topkgenes from h5ad
    '''
    house_keeping_genes = pd.read_csv(housekeeping_genes, index_col = 0)
    top_k_genes = loki.preprocess.generate_gene_df(h5ad, house_keeping_genes)
    return(top_k_genes)


def encode_mtx(top_k_genes, outdir, name, model_path = MODEL_PATH) -> None:
    '''
    for align, anno, decompose
    '''    
    model, preprocess, tokenizer = loki.utils.load_model(model_path, device = DEVICE)
    mtx_embeddings = loki.utils.encode_text_df(model, tokenizer, top_k_genes, 'label', device = DEVICE)
    
    df = pd.DataFrame(mtx_embeddings.cpu().numpy().T)
    df.to_csv(f'{outdir}/{name}_mtx_embeddings.csv', header=False, float_format='%.6f')  # for align, anno
    #np.save(f'{outdir}/{name}_mtx_embeddings.npy', df) 
    torch.save(mtx_embeddings, f'{outdir}/{name}_mtx_embeddings.pt')
    df.columns = top_k_genes.index
    df.to_csv(f'{outdir}/{name}_mtx_embeddings_header.csv', float_format='%.6f')  # for decompose


def encode_image(space_input, coord, outdir, model_path = MODEL_PATH) -> None:
    '''
    for align, anno, decompose
    '''
    ad = sc.read_h5ad(f'{space_input}/valid_spots.h5ad')
    ad.var_names_make_unique()
    print(ad.obs_names[0:3])
    coord = pd.read_csv(coord, index_col=0)
    img = Image.open(f'{space_input}/spatial/tissue_hires_image.png')
    img_array = np.asarray(img)

    # segment
    patch_dir = f'{outdir}/spots_images'
    img_list = os.listdir(patch_dir)
    patch_paths = [os.path.join(patch_dir, fn) for fn in img_list]

    model, preprocess, tokenizer = loki.utils.load_model(model_path, device = DEVICE)
    image_embeddings = loki.utils.encode_images(model, preprocess, patch_paths, device = DEVICE)
    df = pd.DataFrame(image_embeddings.cpu().numpy().T)
    df.columns = [x.replace("_hires.png", "") for x in img_list]
    print(df.columns[0:3])
    df = df.reindex(columns = ad.obs_names)
    print(df.columns[0:3])
    df.to_csv(f'{outdir}/space_image_embeddings.csv', header=False, float_format='%.6f')
    df.to_csv(f'{outdir}/space_image_embeddings_header.csv', float_format='%.6f')
    torch.save(image_embeddings, f'{outdir}/space_image_embeddings.pt')

    ## loki segment
    #patch_dir = f'{outdir}/spots_images_loki'
    #loki.preprocess.segment_patches(img_array, coord, patch_dir)
    #img_list = os.listdir(patch_dir)
    #patch_paths = [os.path.join(patch_dir, fn) for fn in img_list]

    #model, preprocess, tokenizer = loki.utils.load_model(model_path, device = DEVICE)
    #image_embeddings = loki.utils.encode_images(model, preprocess, patch_paths, device = DEVICE)
    #df = pd.DataFrame(image_embeddings.cpu().numpy().T)
    #df.to_csv(f'{outdir}/space_image_embeddings_loki.csv', header=False, float_format='%.6f')
    #torch.save(image_embeddings, f'{outdir}/space_image_embeddings_loki.pt')
    #df.columns = img_list
    #df.to_csv(f'{outdir}/space_image_embeddings_header_loki.csv', float_format='%.6f')


#def encode_sc_text(sc_h5ad, housekeeping_genes, outdir) -> None:
#    '''
#    for decompose
#    '''
#    ad = sc.read_h5ad(sc_h5ad)
#    ad.var_names_make_unique()
#    house_keeping_genes = pd.read_csv(housekeeping_genes, index_col = 0)
#    top_k_genes = loki.preprocess.generate_gene_df(ad, house_keeping_genes)

#    model, preprocess, tokenizer = loki.utils.load_model(MODEL_PATH, device = DEVICE)
#    mtx_embeddings = loki.utils.encode_text_df(model, tokenizer, top_k_genes, 'label', device = DEVICE)
#    df = pd.DataFrame(mtx_embeddings.cpu().numpy().T)
#    df.columns = ad.obs.index
#    df.to_csv(f'{outdir}/sc_mtx_embeddings_header.csv', float_format='%.6f')


def finetune(top_k_genes, spots_images_dir, outdir, spname):
    #ad = sc.read_h5ad(f'{space_input}/valid_spots.h5ad')
    #ad.var_names_make_unique()
    #house_keeping_genes = pd.read_csv(housekeeping_genes, index_col = 0)
    #top_k_genes = loki.preprocess.generate_gene_df(ad, house_keeping_genes)
    top_k_genes['img_idx'] = top_k_genes.index + '_hires'
    top_k_genes['img_path'] = f'{spots_images_dir}/' + top_k_genes.index + '_hires.png'

    top_k_genes.to_csv(f"{outdir}/finetune.csv")

    train_csv = f"{outdir}/finetune.csv"
    name = f'finetune_{spname}'

    train_command = [
        'python', '-m', 'open_clip_train.main',
        '--name', name,
        '--save-frequency', '5',
        '--zeroshot-frequency', '10',
        '--report-to', 'none',   # 'wandb'
        '--train-data', train_csv,
        '--csv-img-key', 'img_path',
        '--csv-caption-key', 'img_idx',  # 'label'
        '--warmup', '10',
        '--batch-size', '64',
        '--lr', '5e-6',
        '--wd', '0.1',
        '--epochs', '10',
        '--workers', '16',
        '--model', MODEL_NAME,
        '--csv-separator', ',',
        '--pretrained', MODEL_PATH,
        '--lock-text-freeze-layer-norm',
        '--lock-image-freeze-bn-stats',
        '--coca-caption-loss-weight','0',
        '--coca-contrastive-loss-weight','1',
        '--val-frequency', '10',
        '--aug-cfg', 'color_jitter=(0.32, 0.32, 0.32, 0.08)', 'color_jitter_prob=0.5', 'gray_scale_prob=0'
    ]
    subprocess.run(train_command)
