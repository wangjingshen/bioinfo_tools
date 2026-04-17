source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/compass/script/plot_update_use_R.R \
    --wilcox_res outdir/wilcox_res.tsv \
    --subsystem Glycolysis/gluconeogenesis \
    --xlab_name KD_vs_WT \
    --outdir outdir \
    --size 2
