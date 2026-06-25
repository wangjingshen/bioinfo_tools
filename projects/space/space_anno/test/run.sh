source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script_dev/space_anno/script/space_anno.R \
    --space_input /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/loki/r5/brain_10X_1/brain_10X_1/space_input/ \
    --sc /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/loki/r5/sc/mus_brain.rds \
    --outdir outdir