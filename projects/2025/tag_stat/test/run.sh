source activate r4.1_env


Rscript /SGRNJ06/randd/USER/wangjingshen/script/tag_stat/script/analysis.R \
    --rds /SGRNJ06/randd/PROJECT/RD20102301_DZH/P24010508_AiBo/20241024_RNA/singleR_annotation/barcodes_9_RNA/barcodes_9_RNA.rds \
    --version singleR \
    --tag_file /SGRNJ06/randd/PROJECT/RD20102301_DZH/P24010508_AiBo/20241024_Target_SNR1.5/barcodes_9_Target/03.count_tag/barcodes_9_Target_tsne_tag.tsv \
    --spname test \
    --outdir outdir