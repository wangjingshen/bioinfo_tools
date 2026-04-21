source activate pyscenic_env

python /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/script/pyscenic.py \
    --input outdir/pyscenic_outdir/pbmc_anno.h5ad \
    --mode h5ad \
    --species human \
    --threads 10 \
    --outdir outdir/pyscenic_outdir \
    #--subset AR,E2F1,SOX4,KLF8 \