source activate jsr4.1

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2022/doublets_stat/script/doublets_stat.R \
    --mix /SGRNJ03/randd/RD20102301_DZH/huashan_HBV/20220207/doublets/H0118_3_HBV_25_R1_ZL.mix_Species.txt \
    --zl_matrix /SGRNJ03/randd/RD20102301_DZH/huashan_HBV/20220207/H0118_3_HBV_25_R1_ZL/05.count/H0118_3_HBV_25_R1_ZL_matrix_10X/ \
    --outdir outdir
