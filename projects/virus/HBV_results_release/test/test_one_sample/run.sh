source activate cellrank

python /SGRNJ06/randd/USER/wangjingshen/script/2023/HBV_results_release/script/analysis.py \
  --rds /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/rds/2022-11-10_HBV-220118_sn_add.rds \
  --prefix HBV_220118_1_EZ_ZL \
  --fj_path /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-12-29_celescope/celescope/FJ/HBV_220118_FJ/ \
  --consensus_path /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-12-30_virus_gene_consensus/gene_consensus/HBV_220118_FJ/ \
  --outdir results \
  --name HBV_220118_FJ