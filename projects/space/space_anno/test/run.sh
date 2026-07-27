source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/space_anno/scripts/pipeline.py \
    --space_dir /SGRNJ06/randd/PROJECT/IR_group/SC_Spatial-transcriptome/0604_Mus_OCT_Nlib/kidney/JC_1/0604_Mus_Kidney_OCT_Nlib/ \
    --sc /SGRNJ07/Standard_Analysis/celelens2local/202606291726_RD24012902_B1/majordataset/RD24012902_B1_major_dataset.rds \
    --score_filter 0 \
    --name Mus_Kidney \
    --outdir outdir