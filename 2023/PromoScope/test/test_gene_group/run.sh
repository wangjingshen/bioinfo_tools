source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/PromoScope/script/gene_group.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/project/P24092403/B1/r38/data/outdir_Main/sweet_tag.rds \
    --cluster Fibroblasts,MuralCells \
    --gene MUC1,MUC2 \
    --outdir outdir