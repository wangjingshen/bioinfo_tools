source activate pyscenic_env

python /SGRNJ06/randd/USER/wangjingshen/script/2023/seurat_scanpy_trans/script/trans.py \
    --rds /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/test/outdir/rds_to_h5ad/pbmc10k.rds \
    --mode rds \
    --outdir outdir \
    --name seurat_scanpy