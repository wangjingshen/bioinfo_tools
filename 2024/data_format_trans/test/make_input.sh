source activate r4.1_env

# seurat -> anndata
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/singleR/test/outdir/seurat_singleR.rds \
    --from_type seurat \
    --to_type anndata \
    --outdir input_dir

# seurat -> sce
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/singleR/test/outdir/seurat_singleR.rds \
    --from_type seurat \
    --to_type sce \
    --outdir input_dir

# seurat -> loom
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/singleR/test/outdir/seurat_singleR.rds \
    --from_type seurat \
    --to_type loom \
    --outdir input_dir
