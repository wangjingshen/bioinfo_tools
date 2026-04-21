# https://stuartlab.org/signac/articles/pbmc_vignette.html
# https://stuartlab.org/signac/articles/mouse_brain_vignette.html

suppressWarnings(suppressMessages({
    library(Signac)
    library(Seurat)
    library(EnsDb.Hsapiens.v75)
    library(EnsDb.Mmusculus.v79)
    library(tidyverse)
    library(SingleR)
    library(argparser)
}))
set.seed(1)
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--analysis_dir", help="analysis dir, such as: /SGRNJ06/randd/PROJECT/scATAC/20230806_mc_959595/")
argv <- add_argument(argv,"--mod", help="analysis mod, doublets or human or mouse")
argv <- add_argument(argv,"--peak_h5", help="peak h5")
argv <- add_argument(argv,"--metadata", help="metadata")
argv <- add_argument(argv,"--fragments", help="fragments")
argv <- add_argument(argv,"--sname", help="sname")
argv <- add_argument(argv,"--resolution", help="resolution, Default: 0.8")
#argv <- add_argument(argv,"--outdir", help="outdir")
argv <- parse_args(argv)

mod <- ifelse(is.na(argv$mod), "doublets", argv$mod)
sname <- argv$sname
if(is.na(argv$peak_h5) & is.na(argv$metadata) & is.na(argv$fragments)){
    peak_h5 <- paste0(argv$analysis_dir, "/", sname, "/03.atac/", sname ,"/outs/filtered_peak_bc_matrix.h5")
    metadata <- paste0(argv$analysis_dir, "/", sname, "/03.atac/", sname ,"/outs/singlecell.csv")
    fragments <- paste0(argv$analysis_dir, "/", sname, "/03.atac/", sname ,"/outs/fragments.tsv.gz")
}else{
    peak_h5 <- argv$peak_h5
    metadata <- argv$metadata
    fragments <- argv$fragments
}
resolution <- ifelse(is.na(argv$resolution), 0.8, as.numeric(argv$resolution))
outdir <- paste0(argv$sname, "_outdir/")
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}


## functions ----
# filter threshold ----
function_clustering <- function(data){
    data <- data %>% 
        RunTFIDF(verbose = F) %>%
        FindTopFeatures(min.cutoff = 'q0', verbose = F) %>%
        RunSVD(verbose = F) %>%
        RunUMAP(reduction = 'lsi', dims = 2:30, verbose = F) %>%
        FindNeighbors(reduction = 'lsi', dims = 2:30, verbose = F) %>%
        FindClusters(algorithm = 3, resolution = 0.8, verbose = FALSE)

    gene.activities <- GeneActivity(data)  # compute gene activities
    data[['RNA']] <- CreateAssayObject(counts = gene.activities)  # add the gene activity matrix to the Seurat object as a new assay
    data <- NormalizeData(object = data, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(data$nCount_RNA))
    return(data)
}

function_singleR <- function(data, species, resolution){
    if(species == "human"){
        ref = readRDS("/SGRNJ03/randd/user/wangjingshen/share/singleR_ref_rds/HumanPrimaryCellAtlasData.rds")
    }
    if(species == "mouse"){
        ref = readRDS("/SGRNJ03/randd/user/wangjingshen/share/singleR_ref_rds/MouseRNAseqData.rds")
    }

    anno <- SingleR(test = data@assays$RNA@data, ref = ref, 
                    clusters = unlist(data@meta.data[[paste0("peaks_snn_res.", resolution)]]),
                    assay.type.test=1, labels = ref$label.main)
    data$cell_type_singleR <- plyr::mapvalues(x = data@meta.data[[paste0("peaks_snn_res.", resolution)]], from = row.names(anno), to = anno$labels)
    return(data)
}

function_peak_cluster <- function(peak_mtx, species, chr_anno){
    if(species == 'human'){
        data <- peak_mtx
        genome <- "hg19"
        # extract gene annotations from EnsDb
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75, verbose = F)
        # change to UCSC style since the data was mapped to hg19
        seqlevels(annotations) <- paste0(chr_anno, seqlevels(annotations))
        genome(annotations) <- "hg19"

        # threshold
        NS_threshold <- 4  # nucleosome_signal
        high_TSS_threshold <- 3
        # filter threshold
        nCount_peaks_min_threshold <- 3000 
        nCount_peaks_max_threshold <- 30000
        pct_reads_in_peaks_threshold <- 15
        blacklist_ratio_threshold <- 0.05
        nucleosome_signal_threshold <- 4
        TSS_enrichment_threshold <- 3
    }

    if(species == 'mouse'){
        data <- peak_mtx
        genome <- "mm10"
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = F)
        seqlevels(annotations) <- paste0(chr_anno, seqlevels(annotations))  ## for doublets
        genome(annotations) <- "mm10"

        # threshold
        NS_threshold <- 4
        high_TSS_threshold <- 2
        # filter threshold
        peak_region_fragments_min_threshold <- 3000
        peak_region_fragments_max_threshold <- 100000
        pct_reads_in_peaks_threshold <- 40
        blacklist_ratio_threshold <- 0.025
        nucleosome_signal_threshold <- 4
        TSS_enrichment_threshold <- 2
    }

    assay <- CreateChromatinAssay(counts = data, sep = c(":", "-"), genome = genome, fragments = fragments, min.cells = 1, verbose = F)
    metadata <- read.csv(file = metadata, header = TRUE, row.names = 1)
    seurat <- CreateSeuratObject(counts = assay, assay = 'peaks', project = 'ATAC', meta.data = metadata)
    Annotation(seurat) <- annotations     # add the gene information to the object
    seurat <- NucleosomeSignal(object = seurat, verbose = F)
    seurat$nucleosome_signal[ is.na(seurat$nucleosome_signal) ] <- 0    # NA to 0
    seurat$nucleosome_signal[ !is.finite(seurat$nucleosome_signal) ] <- 0   # inf to 0
    seurat <- TSSEnrichment(seurat, fast = FALSE, verbose = F)
    seurat$pct_reads_in_peaks <- seurat$peak_region_fragments / seurat$passed_filters * 100
    seurat$blacklist_ratio <- seurat$blacklist_region_fragments / seurat$peak_region_fragments

    # QC plot ----
    # nCount_peaks_TSS_enrichment
    DensityScatter(seurat, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
    ggsave(paste0(outdir, "/", species, "_nCount_peaks_TSS_enrichment.png"), width = 7, height = 5)
    # high tss
    seurat$high.tss <- ifelse(seurat$TSS.enrichment > high_TSS_threshold, 'TSS_high', 'TSS_low')
    TSSPlot(seurat, group.by = 'high.tss') + NoLegend()
    ggsave(paste0(outdir, "/", species, "_high_TSS_enrichment.png"), width = 7, height = 5)
    # nucleosome_signal
    seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > NS_threshold, 'NS_high', 'NS_low')
    #print(table(seurat$nucleosome_group))
    if(length(unique(seurat$nucleosome_group))!=1 & min(table(seurat$nucleosome_group)) > 100){
        FragmentHistogram(object = seurat, group.by = 'nucleosome_group')
        ggsave(paste0(outdir, "/", species, "_nucleosome_group.png"), width = 7, height = 5)
    }else{
        print(table(seurat$nucleosome_group))
        print("Nucleosome_group values low, not plot.")
    }

    # vlnplot
    VlnPlot(object = seurat, features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'), pt.size = 0, ncol = 5)
    ggsave(paste0(outdir, "/", species, "_QC.png"), width = 15, height = 5)

    # not filter ----
    # clustering
    seurat_raw <- seurat
    seurat <- function_clustering(seurat)
    
    DimPlot(object = seurat, label = TRUE)
    ggsave(paste0(outdir, "/", species, "_clusters.png"), width = 7, height = 5)
    seurat <- function_singleR(seurat, species, resolution)
    DimPlot(seurat, reduction = "umap", group.by = "cell_type_singleR", label = T)
    ggsave(paste0(outdir, "/", species, "_cell_type.png"), width = 7, height = 5)
    
    saveRDS(seurat, paste0(outdir, "/", species, ".rds"))

    # filter ----
    if(species == "human"){
        stat_df <- data.frame(condition = c(paste0("nCount_peaks > ", nCount_peaks_min_threshold), 
                                            paste0("nCount_peaks < ", nCount_peaks_max_threshold),
                                            paste0("pct_reads_in_peaks > ", pct_reads_in_peaks_threshold), 
                                            paste0("blacklist_ratio < ", blacklist_ratio_threshold), 
                                            paste0("nucleosome_signal < ", nucleosome_signal_threshold), 
                                            paste0("TSS.enrichment > ", TSS_enrichment_threshold),
                                            "all_condition"),
                              positive_cells = c(sum(seurat_raw$nCount_peaks > nCount_peaks_min_threshold), 
                                                 sum(seurat_raw$nCount_peaks < nCount_peaks_max_threshold),
                                                 sum(seurat_raw$pct_reads_in_peaks > pct_reads_in_peaks_threshold),
                                                 sum(seurat_raw$blacklist_ratio < blacklist_ratio_threshold),
                                                 sum(seurat_raw$nucleosome_signal < nucleosome_signal_threshold), 
                                                 sum(seurat_raw$TSS.enrichment > TSS_enrichment_threshold),
                                                 sum(seurat_raw$nCount_peaks > nCount_peaks_min_threshold & 
                                                     seurat_raw$nCount_peaks < nCount_peaks_max_threshold & 
                                                     seurat_raw$pct_reads_in_peaks > pct_reads_in_peaks_threshold & 
                                                     seurat_raw$blacklist_ratio < blacklist_ratio_threshold &
                                                     seurat_raw$nucleosome_signal < nucleosome_signal_threshold & 
                                                     seurat_raw$TSS.enrichment > TSS_enrichment_threshold)),
                              total_cells = ncol(seurat_raw)) %>%
                              mutate(positive_percent = paste0(round(100*positive_cells/total_cells,2), " %"))
        write.table(stat_df, paste0(outdir, "/", species, "_filter_stat.xls"), sep="\t", row.names = F, quote = F)

        #if(stat_df[stat_df$condition == "all_condition", "positive_cells"] < 50){
        #    print(paste0("Number of cells after filtration: ", stat_df[stat_df$condition == "all_condition", "positive_cells"], " , too low"))
        #    print("Not subset.")
        #}else{
        #    seurat_subset <- subset(seurat_raw, subset = nCount_peaks > nCount_peaks_min_threshold &
        #                                                 nCount_peaks < nCount_peaks_max_threshold &
        #                                                 pct_reads_in_peaks > pct_reads_in_peaks_threshold & 
        #                                                 blacklist_ratio < blacklist_ratio_threshold & 
        #                                                 nucleosome_signal < nucleosome_signal_threshold & 
        #                                                 TSS.enrichment > TSS_enrichment_threshold)
        #    print(paste0("Remove ", ncol(seurat_raw) - ncol(seurat_subset), " cells, remaining ", ncol(seurat_subset), " cells"))
        #    seurat_subset <- function_clustering(seurat_subset)

        #    DimPlot(object = seurat_subset, label = TRUE)
        #    ggsave(paste0(outdir, "/", species, "_filter_clusters.png"), width = 7, height = 5)
        #    seurat_subset <- function_singleR(seurat_subset, species, resolution)
        #    DimPlot(seurat_subset, reduction = "umap", group.by = "cell_type_singleR", label = T)
        #    ggsave(paste0(outdir, "/", species, "_filter_cell_type.png"), width = 7, height = 5)

        #    saveRDS(seurat, paste0(outdir, "/", species, "_filter.rds"))
        #}


    }else{ # "mouse"
        stat_df <- data.frame(condition = c(paste0("peak_region_fragments > ", peak_region_fragments_min_threshold), 
                                            paste0("peak_region_fragments < ", peak_region_fragments_max_threshold),
                                            paste0("pct_reads_in_peaks > ", pct_reads_in_peaks_threshold),
                                            paste0("blacklist_ratio < ", blacklist_ratio_threshold),
                                            paste0("nucleosome_signal < ", nucleosome_signal_threshold),
                                            paste0("TSS.enrichment > ", TSS_enrichment_threshold),
                                            "all_condition"),
                              positive_cells = c(sum(seurat_raw$peak_region_fragments > peak_region_fragments_min_threshold), 
                                                 sum(seurat_raw$peak_region_fragments < peak_region_fragments_max_threshold),
                                                 sum(seurat_raw$pct_reads_in_peaks > pct_reads_in_peaks_threshold),
                                                 sum(seurat_raw$blacklist_ratio < blacklist_ratio_threshold),
                                                 sum(seurat_raw$nucleosome_signal < nucleosome_signal_threshold), 
                                                 sum(seurat_raw$TSS.enrichment > TSS_enrichment_threshold),
                                                 sum(seurat_raw$peak_region_fragments > peak_region_fragments_min_threshold & 
                                                     seurat_raw$peak_region_fragments < peak_region_fragments_max_threshold &
                                                     seurat_raw$pct_reads_in_peaks > pct_reads_in_peaks_threshold & 
                                                     seurat_raw$blacklist_ratio < blacklist_ratio_threshold & 
                                                     seurat_raw$nucleosome_signal < nucleosome_signal_threshold & 
                                                     seurat_raw$TSS.enrichment > TSS_enrichment_threshold)),
                              total_cells = ncol(seurat_raw)) %>%
                              mutate(positive_percent = paste0(round(100*positive_cells/total_cells,2), " %"))
        write.table(stat_df, paste0(outdir, "/", species, "_filter_stat.xls"), sep="\t", row.names = F, quote = F)

        #if(stat_df[stat_df$condition == "all_condition", "positive_cells"] < 50){
        #    print(paste0("Number of cells after filtration: ", stat_df[stat_df$condition == "all_condition", "positive_cells"], ", too low"))
        #    print("Not subset.")
        #}else{
        #    seurat_subset <- subset(seurat_raw, subset = peak_region_fragments > peak_region_fragments_min_threshold & 
        #                                                 peak_region_fragments < peak_region_fragments_max_threshold & 
        #                                                 pct_reads_in_peaks > pct_reads_in_peaks_threshold & 
        #                                                 blacklist_ratio < blacklist_ratio_threshold & 
        #                                                 nucleosome_signal < nucleosome_signal_threshold & 
        #                                                 TSS.enrichment > TSS_enrichment_threshold)
        #    print(paste0("Remove ", ncol(seurat_raw) - ncol(seurat_subset), " cells, remaining ", ncol(seurat_subset), " cells"))
        #    seurat_subset <- function_clustering(seurat_subset)

        #    DimPlot(object = seurat_subset, label = TRUE)
        #    ggsave(paste0(outdir, "/", species, "_filter_clusters.png"), width = 7, height = 5)
        #    seurat_subset <- function_singleR(seurat_subset, species, resolution)
        #    DimPlot(seurat_subset, reduction = "umap", group.by = "cell_type_singleR",label = T)
        #    ggsave(paste0(outdir, "/", species, "_filter_cell_type.png"), width = 7, height = 5)

        #    saveRDS(seurat, paste0(outdir, "/", species, "_filter.rds"))
        #}
    }
}

# single species analysis ----
if(mod == "human"){
    peak_mtx <- Read10X_h5(peak_h5)
    function_peak_cluster(peak_mtx, species="human", chr_anno = 'chr')
}

if(mod == "mouse"){
    peak_mtx <- Read10X_h5(peak_h5)
    function_peak_cluster(peak_mtx, species="mouse", chr_anno = 'chr')
}

## doublets analysis ----
if(mod == "doublets"){
    peak_mtx <- Read10X_h5(peak_h5)
    # split species --
    hs_percent <- colSums(peak_mtx[grep("GRCh38_",row.names(peak_mtx)),])/colSums(peak_mtx)
    mm_percent <- colSums(peak_mtx[grep("mm10_",row.names(peak_mtx)),])/colSums(peak_mtx)
    hs_barcodes <- colnames(peak_mtx)[ hs_percent >= 0.8 & mm_percent <= 0.2] 
    mm_barcodes <- colnames(peak_mtx)[ mm_percent >= 0.8 & hs_percent <= 0.2]
    doublets_barcodes <- colnames(peak_mtx)[ hs_percent < 0.8 & mm_percent < 0.8]

    print(paste0("hs barocdes: ", length(hs_barcodes)))
    print(paste0("mm barcodes: ", length(mm_barcodes)))
    print(paste0("doublets: ", length(doublets_barcodes)))

    function_peak_cluster(peak_mtx[, hs_barcodes], "human", "GRCh38_chr")
    function_peak_cluster(peak_mtx[, mm_barcodes], "mouse", 'mm10_chr')
}
print("Done.")