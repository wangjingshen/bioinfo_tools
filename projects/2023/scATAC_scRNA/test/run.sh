source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/2023/scATAC_scRNA/script/analysis.R \
    --fragments /SGRNJ06/randd/PROJECT/scATAC/summary_10X_data_analysis/mouse/20231206_mouse_embryo/E11_5WT_XA10X_outdir/outs/fragments.tsv.gz,/SGRNJ06/randd/PROJECT/scATAC/self_pipe/mouse_embryo_merge/mouse_embryo/02.atac/Result/Mapping/mouse_embryo/fragments_corrected_dedup_count.tsv.gz \
    --peak_bed /SGRNJ06/randd/PROJECT/scATAC/summary_10X_data_analysis/mouse/20231206_mouse_embryo/E11_5WT_XA10X_outdir/outs/peaks.bed,/SGRNJ06/randd/PROJECT/scATAC/self_pipe/mouse_embryo_merge/mouse_embryo/outs/mouse_embryo_final_peaks.bed \
    --sample pipe_10X,pipe_self \
    --mod 10X,self \
    --species mouse \
    --scRNA /SGRNJ03/randd/lims_result_rd/MultiResult/RD22110301_231130_mus_Embryo/project_batch/2023-12-01/hkcxs2gxo7/father_cluster/batch_1/call-integration/RD22110301_231130_mus_Embryo.diff_PRO.rds \
    --metadata /SGRNJ06/randd/PROJECT/scATAC/summary_10X_data_analysis/mouse/20231206_mouse_embryo/E11_5WT_XA10X_outdir/outs/singlecell.csv,/SGRNJ06/randd/PROJECT/scATAC/self_pipe/mouse_embryo_merge/mouse_embryo/outs/cell_qc_metrics.tsv \
    --filter_bc /SGRNJ06/randd/PROJECT/scATAC/summary_10X_data_analysis/mouse/20231206_mouse_embryo/E11_5WT_XA10X_outdir/outs/filtered_peak_bc_matrix/barcodes.tsv,test \
    --resolution 0.8 \
    --plot_region "chr14-99700000-99760000" \
    --outdir outdir \
#    --plot_cluster Neurons \
#    --plot_cluster_region "chr11-98279644-98379645" \    