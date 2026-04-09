source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/script/seuratGtf/script/pipeline.py \
    --matrix_10X /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/ffpe/r2/rna/LFPS260302PXL/outs/filtered/,/SGRNJ06/randd/USER/wangjingshen/rd_project/2026/ffpe/r2/rna/LFPS260302YWL/outs/filtered/ \
    --spname Subcutaneous_tumor,In_situ_tumor \
    --gname Subcutaneous_tumor,In_situ_tumor \
    --species mouse \
    --resolution 0.8 \
    --gtf /SGRNJ06/randd/public/genome/rna/mmu/mmu_ensembl_110_nofilter/Mus_musculus.GRCm39.110.gtf \
    --gene_type protein_coding \
    --outdir outdir \
