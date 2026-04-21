source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/QC/script/qc_vlnplot.R \
    --rds /SGRNJ06/randd/PROJECT/RD20102301_DZH/PK22011902_huashanHBV/20220909/HS0826-1/06.analysis/HS0826-1.rds \
    --species human \
    --sname HS0826 \
    --outdir outdir
