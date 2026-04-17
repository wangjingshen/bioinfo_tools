source activate scenic

python /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/script/pyscenic_regulon_plot.py \
    --pyscenic_metadata_loom pyscenic_plot/pyscenic_output_add_metadata.loom \
    --h5ad pyscenic_outdir/pbmc_anno.h5ad \
    --network pyscenic_outdir/cytoscape_network.txt \
    --node pyscenic_outdir/cytoscape_node.txt \
    --regulons BCL11A,EOMES \
    --ncol 2 \
    --outdir outdir/pyscenic_plot
