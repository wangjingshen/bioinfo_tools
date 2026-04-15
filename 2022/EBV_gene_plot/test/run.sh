source activate jsr4.1

Rscript /SGRNJ06/randd/USER/wangjingshen/script/EBV_gene_plot/EBV_gene_plot.R \
    --zl_matrix /SGRNJ03/randd/RD20102301_DZH/virus_panel/20211110_6/E1026_Rji_3T3_SN424211_LGR-3ZL/06.count/E1026_Rji_3T3_SN424211_LGR-3ZL_matrix_10X/ \
    --virus_matrix /SGRNJ06/randd/PROJECT/RD20102301_DZH/data_for_display/EBV_panel/celescope1.13.0/20221110/SQ_Raji_3T3_SN424211_LGR-3TS/09.count/virus_matrix/ \
    --outdir .
