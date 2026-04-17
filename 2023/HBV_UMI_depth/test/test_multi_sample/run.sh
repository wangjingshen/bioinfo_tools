source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/HBV_UMI_depth/script/plot.R \
    --virus_id_file ../get_virus_bam/R220121031/R220121031_virus_id.tsv,../get_virus_bam/R220121032/R220121032_virus_id.tsv,../get_virus_bam/R220124006/R220124006_virus_id.tsv,../get_virus_bam/R220124007/R220124007_virus_id.tsv \
    --sample_name R220121031,R220121032,R220124006,R220124007 \
    --total_reads 41435652,23128536,6147385,6313889 \
    --downsample_mode 10