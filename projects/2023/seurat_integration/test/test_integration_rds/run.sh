source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/seurat_integration/script/integration.R \
    --mode rds \
    --rds /SGRNJ06/randd/USER/wangjingshen/script/2023/seurat_integration/data_local/pancreas.rds \
    --batch_var tech \
    --celltype_var celltype