source activate compass

python ../script/compass_plot.py \
    --reaction_penalties compass_outdir/reactions.tsv \
    --metadata ../data_local/reaction_metadata.csv \
    --group_info input/group_info.tsv \
    --subsystem Glycolysis/gluconeogenesis \
    --outdir compass_plot_outdir \
    --name Glycolysis_gluconeogenesis
