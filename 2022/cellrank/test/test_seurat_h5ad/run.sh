python /SGRNJ03/randd/user/wangjingshen/project/cellrank/cellrank_pipeline/run_cellrank.py \
  --sc_file /SGRNJ03/randd/user/wangjingshen/project/cellrank/10X_1k_mouse_brain/data/neuron_seurat.h5ad \
  --loom_file /SGRNJ03/randd/user/wangjingshen/project/cellrank/10X_1k_mouse_brain/data/neuron_1k_v3_possorted_genome_sorted_UY8NC.loom \
  --cluster_key cell_type \
  --n_terminal_states 2 \
  --n_initial_states 1 \
  --outdir results/ \
  --save_h5ad T
  #--driver_genes driver_genes.csv \
