source activate scenic

python /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/script/h5ad_anno.py \
    --h5ad /SGRNJ06/randd/USER/wangjingshen/script/pyscenic/test/test_pbmc/pyscenic_outdir/pbmc.h5ad \
    --anno input/anno.csv \
    --mode csv \
    --outdir pyscenic_outdir \
    --name pbmc_anno
