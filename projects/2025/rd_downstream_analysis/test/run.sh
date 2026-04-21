source activate /SGRNJ06/randd/USER/wangjingshen/soft/miniforge3/envs/CeleScope2.6.0_probe

python /SGRNJ06/randd/USER/wangjingshen/script/rna_down/script/pipeline.py \
    --mapfile ./mapfile \
    --rna_mapfile /SGRNJ06/randd/PROJECT/RD20102301_DZH/dzh_test/cuishuang/Paris_Brain/20250630_RNA.mapfile \
    --genomeDir /SGRNJ06/randd/public/genome/rna/celescope_v2/hs \
    --species human \
    --probe_file /SGRNJ06/randd/PROJECT/RD20102301_DZH/dzh_test/cuishuang/Paris_Brain/Paris_Brain_probe.fasta \
    --probe_file_mode JMML
