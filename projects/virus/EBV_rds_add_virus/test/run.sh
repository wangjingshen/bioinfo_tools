source activate jsr4.1

Rscript ../script/EBV_rds_add_virus/analysis.R \
  --rds /SGRNJ03/randd/Integratedanalysis/display/panel/EBV/RAJI3T3/20211123/20220104/amulti/RAJI3T3.diff_PRO.rds \
  --virus_mtx /SGRNJ06/randd/PROJECT/RD20102301_DZH/EBV_panel/outward_display_data/20221220/SQ_Raji_3T3_SN424211_LGR-3TS/09.count/SQ_Raji_3T3_SN424211_LGR-3TS_virus_matrix/ \
  --prefix RAJI3T3 \
  --name RAJI3T3_virus