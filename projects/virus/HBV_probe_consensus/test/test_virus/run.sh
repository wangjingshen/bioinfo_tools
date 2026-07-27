source activate cellrank

python /SGRNJ06/randd/USER/wangjingshen/script/HBV_probe_consensus/HBV_probe_virus_bam_split.py \
    --fj_path /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-07-13_zh_probe/R211209011_fj_celescope/H1207_2_HBV_AD38_FJ/ \
    --probe_anno /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-07-13_zh_probe/R211209011/H1207_2_HBV_AD38_FJ_annoate_result_mis_2.txt \
    --otsu_min_support_read 6
