source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/compass/script/get_input.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/project/datasets/pancreatic_cancer/data/seurat.rds \
    --subset Epithelial_cells \
    --subset_var cell_type_singleR \
    --outdir outdir \
    --name Epithelial_cells