source activate scenic

python /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/script/pyscenic_plot.py \
    --pyscenic_loom pyscenic_outdir/pyscenic_output.loom \
    --h5ad pyscenic_outdir/pbmc_anno.h5ad \
    --threads 10 \
    --ncol 3 \
    --outdir outdir/pyscenic_plot