source activate kb_env

python /SGRNJ06/randd/USER/wangjingshen/script/fasta_to_gtf/script/fasta_to_gtf.py \
    run \
    --genome /SGRNJ06/randd/PROJECT/RD20102301_DZH/dzh_test/cuishuang/P26020701_XTL_GOTN_EBV/XTL_GOTN_EBV_genome \
    --fasta /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/ebv/r2/XTL_EBV_gene.fa \
    --star /SGRNJ/Public/Software/conda_env/celescope2.1.0/bin/STAR-avx2 \
    --prefix XTL_EBV

# blast update EBNA2 
cp /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/ebv/r2/XTL_EBV.bed XTL_EBV_update.bed

python /SGRNJ06/randd/USER/wangjingshen/script/fasta_to_gtf/script/fasta_to_gtf.py \
    update \
    --bed XTL_EBV_update.bed \
    --prefix XTL_EBV