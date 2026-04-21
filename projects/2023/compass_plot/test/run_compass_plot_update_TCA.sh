source activate r4.1_env

Rscript ../script/compass_plot_update.R \
    --wilcox_res compass_plot_outdir/wilcox_res.tsv \
    --subsystem "Citric acid cycle" \
    --xlab_name Th17p_vs_Th17n \
    --annotation_text_size 2 \
    --outdir compass_plot_outdir \
    --name Citric_acid_cycle
