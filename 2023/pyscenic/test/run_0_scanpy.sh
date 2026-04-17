source activate pyscenic_env

python /SGRNJ06/randd/USER/wangjingshen/script/2023/pyscenic/script/scanpy_run.py \
    --mtx /SGRNJ06/randd/USER/wangjingshen/rd_project/pyscenic/data/2023-01-30_test/filtered_feature_bc_matrix/ \
    --species human \
    --threads 10 \
    --n_neighbors 15 \
    --n_pcs 20 \
    --resolution 0.4 \
    --outdir pyscenic_outdir \
    --name pbmc
