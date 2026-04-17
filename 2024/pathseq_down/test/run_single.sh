source activate r4.1_env

# 11577100
Rscript /SGRNJ06/randd/USER/wangjingshen/script/2024/pathseq_down/script/analysis.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/project/P24041601_WuHanKouQiang/B1/r2_analysis/rds/seurat.rds \
    --df_genus /SGRNJ06/randd/USER/wangjingshen/project/P24041601_WuHanKouQiang/B1/r2/N/N_11577100_Target/04.analysis/N_11577100_Target_genus_umi.tsv,/SGRNJ06/randd/USER/wangjingshen/project/P24041601_WuHanKouQiang/B1/r2/T/T_11577100_Target/04.analysis/T_11577100_Target_genus_umi.tsv \
    --subset N_11577100_RNA,T_11577100_RNA \
    --rna_spname N_11577100_RNA,T_11577100_RNA \
    --spname N_11577100_Target,T_11577100_Target \
    --barplot_topn 10 \
    --outdir outdir_11577100 \
    --saveRDS T
