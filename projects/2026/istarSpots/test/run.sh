source activate istar_env

python /SGRNJ06/randd/USER/wangjingshen/script_dev/istarSpots/script/pipeline.py \
    --istar_labels /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/istar/r5_raw/YWL_NJU_ST_Lib/outs/clusters-gene/labels.pickle \
    --dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20260318_tumor/YWL_NJU_ST_Lib/ \
    --spname YWL_NJU_ST_Lib \
    --k 5 \
    --distance_thresh 200 \