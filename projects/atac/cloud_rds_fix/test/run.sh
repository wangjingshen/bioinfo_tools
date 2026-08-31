source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/atac/cloud_rds_fix/scripts/cloud_rds_fix.R \
    --rds /SGRNJ06/randd/PROJECT/RD20102301_DZH/dzh_test/dingixuheng/atac_shouhou/20250305_rna/data_load/cloud_annotation_result.rds \
    --set_group_by_sample F \
    --outdir outdir