source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/PromoScope/script/add_sweet_tag.R \
    --rds /SGRNJ06/randd/USER/wangjingshen/project/P23031708/data/B4/r3/result/Modules/1.ClustAnno/CellLabel/Main/P23031708.diff_PRO.rds \
    --tag_file /SGRNJ06/bioinfo/PROJ04/PROJ_22.Other/P23031708_sweetseq/20230426/NPBMC202304_c3/NPBMC202304/04.count_tag/NPBMC202304_tsne_tag.tsv,/SGRNJ06/bioinfo/PROJ04/PROJ_22.Other/P23031708_sweetseq/20230625/T2DM_c3/T2DM/04.count_tag/T2DM_tsne_tag.tsv,/SGRNJ06/bioinfo/PROJ04/PROJ_22.Other/P23031708_sweetseq/20230703/DKD_c3/DKD/04.count_tag/DKD_tsne_tag.tsv \
    --sname NPBMC,T2DM,DKD \
    --subcluster_rds /SGRNJ06/randd/USER/wangjingshen/project/P23031708/data/B4/r3/result/Modules/1.ClustAnno/CellLabel/TandNK_2/P23031708.diff_PRO.rds \
    --outdir sweet_tag \
    --name outdir