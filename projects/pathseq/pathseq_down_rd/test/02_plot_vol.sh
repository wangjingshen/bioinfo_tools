source /SGRNJ/Public/Software/anaconda3/bin/activate Seuratv3

Rscript /Public/Script/shouhou/Integration_Platform/Function_groupDiff/plot_diff_Volcano.R \
    --diff /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/test/r3/outdir/02.diff_rna/Fusobacterium/degs_Fibroblasts_Fusobacterium_positive_vs_negative.tsv \
    --method pvalue_adj \
    --label F \
    --title positive_vs_negative \
    --logfc 0.25 \
    --pvalue 0.05 \
    --top 10 \
    --outdir outdir/02.diff_rna/Volcano \
    --prefix test \
    --version v1 \
