source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/dgc_tag_plot/script/analysis.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/dgc_seq/r1/result.rds \
    --spname Control,TAG1217,TAG1222 \
    --tsne_tag /SGRNJ06/randd/PROJECT/DGC_seq/20260104_TAG/TAG0523_TAG_1224/03.count_tag/TAG0523_TAG_1224_tsne_tag.tsv,/SGRNJ06/randd/PROJECT/DGC_seq/20260104_TAG/TAG1217_TAG_1224/03.count_tag/TAG1217_TAG_1224_tsne_tag.tsv,/SGRNJ06/randd/PROJECT/DGC_seq/20260104_TAG/TAG1222_TAG_1224/03.count_tag/TAG1222_TAG_1224_tsne_tag.tsv \
    --tag_name New_DGC_barcode \
    --VlnPlot T \
    --FeaturePlot F \
    --outdir outdir
