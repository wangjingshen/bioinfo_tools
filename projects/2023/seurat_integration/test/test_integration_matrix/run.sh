source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/seurat_integration/script/integration.R \
    --mode matrix \
    --matrix_10X /SGRNJ06/randd/USER/wangjingshen/rd_project/pyscenic/data/2023-01-30_test/filtered_feature_bc_matrix/,/SGRNJ06/randd/PROJECT/RD20102301_DZH/P22021901_Becowei_EBV/20220613/SPB0531J/05.count/SPB0531J_filtered_feature_bc_matrix/ \
    --batch_name 10X,singleron \