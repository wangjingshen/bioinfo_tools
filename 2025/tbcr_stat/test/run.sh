source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/tbcr_stat/script/analysis.R \
    --rds /SGRNJ06/randd/PROJECT/RD20102301_DZH/P23062001_ZhongshanEBV/20250120_RNA/singleR_annotation/NPC_BM_RNA/NPC_BM_RNA.rds \
    --mode TCR \
    --anno_mode rd_sr \
    --fj_path /SGRNJ06/randd/PROJECT/RD20102301_DZH/P23062001_ZhongshanEBV/20250121_TCR/ \
    --fj_spname NPC_BM_TCR \
    --rna_spname NPC_BM_RNA \
    --outdir outdir/NPC_BM_TCR \


Rscript /SGRNJ06/randd/USER/wangjingshen/script/tbcr_stat/script/analysis.R \
    --rds /SGRNJ06/randd/PROJECT/RD20102301_DZH/P23062001_ZhongshanEBV/20250120_RNA/singleR_annotation/NPC_BM_RNA/NPC_BM_RNA.rds \
    --mode BCR \
    --anno_mode rd_sr \
    --fj_path /SGRNJ06/randd/PROJECT/RD20102301_DZH/P23062001_ZhongshanEBV/20250121_BCR/ \
    --fj_spname NPC_BM_BCR \
    --rna_spname NPC_BM_RNA \
    --outdir outdir/NPC_BM_BCR \