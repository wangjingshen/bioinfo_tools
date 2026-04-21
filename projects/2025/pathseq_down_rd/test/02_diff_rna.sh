source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/script/diff_rna.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/test/r3/outdir/00.data/data_seurat.rds \
    --genus_analysis Fusobacterium \
    --outdir outdir

Rscript /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/script/diff_rna.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/test/r3/outdir/00.data/data_seurat.rds \
    --genus_analysis total_genus \
    --outdir outdir