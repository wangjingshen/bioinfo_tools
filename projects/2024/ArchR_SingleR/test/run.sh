source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2024/ArchR_SingleR/script/analysis.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/project/10X_ATAC/data/2024-02-21_lung_T/run_pipe/01.Main.anno/0.rds/projHeme3.anno.rds \
    --species human \



Rscript /SGRNJ06/randd/USER/wangjingshen/script/2024/ArchR_SingleR/script/analysis.R \
    --fragments /SGRNJ06/randd/PROJECT/scATAC/20240206_CCRF_A549_828_sc/A0124_2_CCRF_A549_6th_D02SDS_T7_EDTA5030/03.atac/A0124_2_CCRF_A549_6th_D02SDS_T7_EDTA5030/outs/fragments.tsv.gz,/SGRNJ06/randd/PROJECT/scATAC/20240206_CCRF_A549_828_sc/A0124_3_CCRF_A549_7th_D02SDS_T7_EDTA5030/03.atac/A0124_3_CCRF_A549_7th_D02SDS_T7_EDTA5030/outs/fragments.tsv.gz \
    --sname A0124_2_CCRF_A549_6th_D02SDS_T7_EDTA5030,A0124_3_CCRF_A549_7th_D02SDS_T7_EDTA5030 \
    --species human \
    
