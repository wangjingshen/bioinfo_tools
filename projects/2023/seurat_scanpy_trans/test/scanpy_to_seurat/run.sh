source activate pyscenic_env

python /SGRNJ06/randd/USER/wangjingshen/script/2023/seurat_scanpy_trans/script/trans.py \
    --h5ad /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/test/outdir/pyscenic_outdir/pbmc_anno.h5ad \
    --mode h5ad \
    --outdir outdir \
    --name scanpy_seurat