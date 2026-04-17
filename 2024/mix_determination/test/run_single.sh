source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2024/mix_determination/script/mix_determination.R \
    --matrix_10X /SGRNJ06/randd/PROJECT/DGC_seq/20240829/Cellmix_0823/outs/filtered/ \
    --matrix_10X_gene /SGRNJ06/randd/PROJECT/DGC_seq/20240829/Cellmix_0823/outs/filtered/features.tsv.gz \
    --outdir outdir \
    --name test \