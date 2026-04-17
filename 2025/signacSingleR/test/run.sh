source activate r4.1_env
Rscript /SGRNJ06/randd/USER/wangjingshen/script/Signac_singleR/script/analysis.R\
  --atac_dir /SGRNJ06/randd/PROJECT/scATAC/project/self_pipe_new_20250805/human/P25041807_xiamen_deng_human/20250819/LH7C_23269984_X2/,/SGRNJ06/randd/PROJECT/scATAC/project/self_pipe_new_20250805/human/P25041807_xiamen_deng_human/20250819/LH7C_23269984/,/SGRNJ06/randd/PROJECT/scATAC/project/self_pipe_new_20250805/human/P25041807_xiamen_deng_human/20250819/LH7C_23269984_X3/ \
  --sample LH7C_23269984_X2,LH7C_23269984,LH7C_23269984_X3 \
  --group LH7C_23269984_X2,LH7C_23269984,LH7C_23269984_X3 \
  --mod self\
  --rm_batch yes \
  --outdir outdir\
  --species human\
  --resolution 0.8 \
  --saveRDS True
