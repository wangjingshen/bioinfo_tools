source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/HBV_UMI_depth/script/plot.R \
    --virus_id_file ../get_virus_bam/R220124006/R220124006_virus_id.tsv \
    --sample_name R220124006 \
    --total_reads 6147385 \
    --downsample_mode 10