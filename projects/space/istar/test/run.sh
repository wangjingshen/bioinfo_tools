source /SGRNJ01/Public/Software/anaconda3/bin/activate /SGRNJ06/randd/USER/wangjingshen/soft/miniforge3/envs/istar_env/

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/istar/scripts/pipeline.py \
    --dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20251127/Mint_FFPE_96_96_1119/ \
    --image /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20251127/Mint_FFPE_96_96_1119/outs/spatial/figure/20251118_chang.jpeg \
    --spname Mint_FFPE_96_96_1119 \
    --swap_pos T \
    --foreground_method istar \
    --foreground_cluster_method max \
    --cluster_method km \
    --n_cluster 10 \
    --k 3 \
    --distance_thresh 200 \
    --clip \
    --step input,preprocess,getMask,impute,cluster,downstream,output,istarSpots,orgdir
    # --hard_radius    #help='set hard radius'


python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/space/istar/scripts/pipeline.py \
    --dir /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/istar/r6/celescope/PXL_NJU_ST_Lib/ \
    --image /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/istar/r6/celescope/PXL_NJU_ST_Lib/outs/spatial/figure/Subcutaneous_tumor_0.05.jpg \
    --spname PXL_NJU_ST_Lib \
    --swap_pos T \
    --foreground_method in_tissue \
    --foreground_cluster_method max \
    --cluster_method km \
    --n_cluster 10 \
    --k 3 \
    --distance_thresh 200 \
    --clip \
    --step input,preprocess,getMask,impute,cluster,downstream,output,istarSpots,orgdir
    # --hard_radius    #help='set hard radius'

