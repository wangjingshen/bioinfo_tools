source activate r4.1_env

# anndata -> loom
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/test/input_dir/seurat_to_anndata.h5ad \
    --from_type anndata \
    --to_type loom \
    --outdir outdir

# seurat -> loom
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/singleR/test/outdir/seurat_singleR.rds \
    --from_type seurat \
    --to_type loom \
    --outdir outdir

# anndata -> sce
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/test/input_dir/seurat_to_anndata.h5ad \
    --from_type anndata \
    --to_type sce \
    --outdir outdir

# sce -> seurat
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/test/input_dir/seurat_to_sce.rds \
    --from_type sce \
    --to_type seurat \
    --outdir outdir

# loom -> seurat
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/test/input_dir/seurat_to_loom.loom \
    --from_type loom \
    --to_type seurat \
    --outdir outdir

# loom -> cds
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/test/input_dir/seurat_to_loom.loom \
    --from_type loom \
    --to_type cds \
    --outdir outdir

# sce -> cds
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/test/input_dir/seurat_to_sce.rds \
    --from_type sce \
    --to_type cds \
    --outdir outdir

# seurat -> cds
Rscript /SGRNJ06/randd/USER/wangjingshen/script/data_format_trans/script/data_format_trans.R \
    --input /SGRNJ06/randd/USER/wangjingshen/script/singleR/test/outdir/seurat_singleR.rds \
    --from_type seurat \
    --to_type cds \
    --outdir outdir