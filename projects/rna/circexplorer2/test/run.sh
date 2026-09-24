source activate celescope3.0.0

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/rna/circexplorer2/scripts/pipeline.py \
    --celescope_dir mouse-testicle/ \
    --name mouse-testicle \
    --refFlat /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/circexplorer2/genome/Mus_musculus.GRCm39.110.refFlat_11col \
    --reference_fa /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/circexplorer2/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --match_window 10
