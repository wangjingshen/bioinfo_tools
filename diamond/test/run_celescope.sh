source activate celescope2.0.7
multi_rna\
  --mapfile mapfile_celescope \
  --STAR_param "--soloStrand Reverse --outReadsUnmapped Fastx" \
  --chemistry customized \
  --pattern C6U12 \
  --whitelist /SGRNJ06/randd/PROJECT/smartseq2-hqq/20250523_fixed/all_hexamers.txt \
  --genomeDir /SGRNJ06/randd/public/genome/rna/celescope_v2/mmu/\
