source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/script/2024/pathseq_down/script/pipeline.py \
    --matrix_10X /SGRNJ06/randd/PROJECT/RD20102301_DZH/P24112214_HJJY_16S/20241220_RNA/WT_3_RNA/outs/filtered/,/SGRNJ06/randd/PROJECT/RD20102301_DZH/P24112214_HJJY_16S/20241220_RNA/lbj_3_RNA/outs/filtered/ \
    --rna_spname WT_3_RNA,lbj_3_RNA \
    --gname WT_3_RNA,lbj_3_RNA \
    --fj_path /SGRNJ06/randd/USER/wangjingshen/project/P24112214_HJJY_16S/B1/r1/report/ \
    --fj_spname WT_3_Target,lbj_3_Target \