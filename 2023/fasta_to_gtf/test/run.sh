source activate /SGRNJ06/randd/USER/wangjingshen/soft/miniforge3/envs/kb_env/

python ../script/fasta_to_gtf.py \
    --genome /SGRNJ06/randd/PROJECT/RD20102301_DZH/P23062001_ZhongshanEBV/EBV_genome \
    --fasta /SGRNJ06/randd/PROJECT/RD20102301_DZH/P23062001_ZhongshanEBV/250303new_18region_gtf/EBV250303new_18region.fasta \
    --prefix EBV250303new_18region
