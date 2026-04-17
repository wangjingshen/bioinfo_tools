source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/script/auto_assign/script/pipeline.py \
    --rds /SGRNJ06/randd/USER/wangjingshen/script_dev/seuratGtf/test/outdir/seurat.rds \
    --spname test \
    --anno_var seurat_clusters \
    --marker_ref /SGRNJ06/randd/USER/wangjingshen/script/auto_assign/marker/mus/bladder.tsv \
    
