source activate cell2loc_env2

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/cell2location/scripts/pipeline.py \
     --sc_ref /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/cell2location/test/test/tutorial/data/sc/sc.h5ad \
     --space_dir /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/cell2location/test/test/tutorial/data/space/ \
     --spname Human_Lymph_Node \
     --sc_max_epochs 250 \
     --N_cells_per_location 30 \
     --detection_alpha 20 \
     --space_max_epochs 30000 \
     --labels_key Subset \
     --batch_key Sample \
     --categorical_covariate_keys 'Method' \
     --cell_count_cutoff 5 \
     --cell_percentage_cutoff2 0.03 \
     --nonz_mean_cutoff 1.12 \
     --step sc_train,space_train,plot
