source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/space_tag/scripts/analysis.R \
    --space_input /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/space_tag/r1/MusBladder_1_OCT_Nlib/ \
    --tag_umi /SGRNJ06/randd/PROJECT/IR_group/SC_Spatial-transcriptome/0709_MusBladder_ND_OCT_Nlib/PBS/DGC/0709_MusBladder_1_TAG/03.count_tag/0709_MusBladder_1_TAG_umi_tag.tsv \
    --outdir outdir