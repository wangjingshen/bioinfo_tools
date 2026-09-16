source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/space_virus/scripts/space_analysis_virus.py \
    --space_dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20260724_NJU_bla/Mus_bla0717_BCG_1_Lib/ \
    --virus_df Mus_bla0717_BCG_1_Lib/06.filter_virus/Mus_bla0717_BCG_1_Lib_filtered_UMI.csv \
    --name Mus_bla0717_BCG_1_Lib

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/space_virus/scripts/space_analysis_virus.py \
    --space_dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20260724_NJU_bla/Mus_bla0717_ENC_1_Lib/ \
    --virus_df Mus_bla0717_ENC_1_Lib/06.filter_virus/Mus_bla0717_ENC_1_Lib_filtered_UMI.csv \
    --name Mus_bla0717_ENC_1_Lib

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/space_virus/scripts/space_analysis_virus.py \
    --space_dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20260724_NJU_bla/Mus_bla0717_PBS_1_Lib/ \
    --virus_df Mus_bla0717_PBS_1_Lib/06.filter_virus/Mus_bla0717_PBS_1_Lib_filtered_UMI.csv \
    --name Mus_bla0717_PBS_1_Lib

    