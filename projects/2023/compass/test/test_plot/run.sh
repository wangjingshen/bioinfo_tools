source activate compass

python /SGRNJ06/randd/USER/wangjingshen/script/2023/compass/script/plot.py \
    --reaction_penalties /SGRNJ06/randd/USER/wangjingshen/project/compass/test/test_scFEA_data/reactions.tsv \
    --metadata /SGRNJ06/randd/USER/wangjingshen/script/compass/data_local/reaction_metadata.csv \
    --group_metadata /SGRNJ06/randd/USER/wangjingshen/script/compass/test/test_plot/input/group_metadata.tsv \
    --subsystem Glycolysis/gluconeogenesis \
    --outdir outdir
