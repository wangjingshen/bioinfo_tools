source /SGRNJ/Public/Software/anaconda3/bin/activate clusterprofiler_update

Rscript /Public/Script/shouhou/Integration_Platform/Function_SigEnrich/clusterprofiler.R \
    --diff /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/test/r3/outdir/02.diff_rna/Fusobacterium/degs_Fibroblasts_Fusobacterium_positive_vs_negative.tsv \
    --species homo_sapiens \
    --updown updown \
    --KEGG TRUE \
    --GO TRUE \
    --Disease FALSE \
    --Reactome F \
    --Wiki FALSE \
    --col 1,2 \
    --FC_cutoff 0.25 \
    --P_cutoff 0.05 \
    --count 3 \
    --cutoff p.adjust \
    --topshow 10 \
    --prefix test \
    --outdir /SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/test/r3/outdir/03.enrich/ \

