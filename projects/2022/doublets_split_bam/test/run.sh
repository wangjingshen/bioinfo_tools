source activate cellrank

python /SGRNJ06/randd/USER/wangjingshen/script/2022/doublets_split_bam/script/doublets_split_bam.py \
    --input_bam /SGRNJ06/randd/USER/wangjingshen/script/doublets_split_bam/test/test_data/H0118_4_polyT_25_R1_FJ/04.star_virus/H0118_4_polyT_25_R1_FJ_virus_Aligned.sortedByCoord.out.bam \
    --doublets_species /SGRNJ03/randd/RD20102301_DZH/huashan_HBV/20220207/doublets/H0118_4_polyT_25_R1_ZL.mix_Species.txt \
    --fj_path /SGRNJ06/randd/USER/wangjingshen/script/doublets_split_bam/test/test_data/H0118_4_polyT_25_R1_FJ/ \
    --otsu_min_support_read 3 \
    --outdir outdir/