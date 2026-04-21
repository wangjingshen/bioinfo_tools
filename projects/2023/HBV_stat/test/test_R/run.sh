source activate jsr4.1

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/HBV_stat/HBV_stat.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/rds/2022-04-02_HBV_220118_1_EZ_2_NEB_ZL.rds \
    --mod same \
    --version v1 \
    --fj_path /SGRNJ03/randd/RD20102301_DZH/huashan_HBV/20220407 \
    --fj_name HBV_220118_1_EZ_FJ,HBV_220118_2_NEB_FJ \
    --prefix HBV_220118_1_EZ_ZL,HBV_220118_2_NEB_ZL
