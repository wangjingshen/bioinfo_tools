source activate loki_env

python /SGRNJ06/randd/USER/wangjingshen/script/loki/script/preprocess.py \
    --dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20251127/Mint_FFPE_96_96_1119/ \
    --spname Mint_FFPE_96_96_1119 \
    --hk_genes /SGRNJ06/randd/USER/wangjingshen/script/loki/data/housekeeping_genes_mus.csv \
    --sc /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/loki/r2/decompose/input/musculus_intestine_subset.h5ad \