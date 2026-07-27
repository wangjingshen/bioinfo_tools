source activate cellrank

python /SGRNJ06/randd/USER/wangjingshen/script/2023/HBV_stat/script/HBV_stat.py \
    --h5ad /SGRNJ06/randd/PROJECT/RD20102301_DZH/HBV_panel/202230428_H0418_rna/H0418_1_Mix1_5M_ZL/06.analysis/H0418_1_Mix1_5M_ZL.h5ad \
    --virus_UMI /SGRNJ06/randd/PROJECT/RD20102301_DZH/HBV_panel/202230428_H0418_virus/H0418_1_Mix1_5M_FJ/07.analysis_virus/H0418_1_Mix1_5M_FJ_UMI_tsne.csv \
    --outdir outdir
