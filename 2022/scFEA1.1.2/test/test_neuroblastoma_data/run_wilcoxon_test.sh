source activate compass

python ../../script/wilcoxon_test.py \
  --flux 02.scFEA_predication/all_module_flux.csv \
  --metadata nb_metadata.tsv \
  --gname Tumour,Endothelium  \
  --supermodule 'Glycolysis' \
  --outdir wilcoxon_test_outdir 
