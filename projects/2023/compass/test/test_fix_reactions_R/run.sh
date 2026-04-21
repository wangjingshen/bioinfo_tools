source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/compass/script/fix_reaction.R \
    --reactions /SGRNJ06/randd/USER/wangjingshen/project/compass/test/test_lad_lsq/reactions.tsv \
    --outdir outdir \
    --name test