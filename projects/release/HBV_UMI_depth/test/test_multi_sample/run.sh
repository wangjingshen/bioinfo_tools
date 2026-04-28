source activate r4.1_env

Rscript ../../script/plot.R \
    --virus_id_file ../get_virus_bam/R220121031/R220121031_virus_id.tsv,../get_virus_bam/R220124006/R220124006_virus_id.tsv \
    --sample_name R220121031,R220124006 \
    --total_reads 41435652,6147385 \
    --downsample_mode 10
