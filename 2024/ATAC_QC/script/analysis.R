# https://stuartlab.org/signac/articles/pbmc_vignette.html
# https://stuartlab.org/signac/articles/mouse_brain_vignette.html

suppressWarnings(suppressMessages({
    library(Signac)
    library(Seurat)
    library(harmony)
    library(EnsDb.Hsapiens.v75)
    library(EnsDb.Mmusculus.v79)
    library(tidyverse)
    library(SingleR)
    library(argparser)
}))
set.seed(1)
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--peak_h5", help="peak h5")
argv <- add_argument(argv,"--metadata", help="metadata")
argv <- add_argument(argv,"--fragments", help="fragments")
argv <- add_argument(argv,"--species", help="species")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--name", help="name")
argv <- parse_args(argv)

peak_h5 <- argv$peak_h5
metadata <- argv$metadata
fragments <- argv$fragments
species <- argv$species
outdir <- argv$outdir
name <- argv$name

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# annotation
if(species == "mouse"){
    genome <- "mm10"
    annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = F))
    seqlevels(annotations) <- paste0("chr", seqlevels(annotations)) 
    genome(annotations) <- "mm10"
}
if(species == "human"){
    genome <- "hg19"
    annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75, verbose = F))
    seqlevels(annotations) <- paste0("chr", seqlevels(annotations)) 
    genome(annotations) <- "hg19"
}

# analysis
data <- Read10X_h5(peak_h5)
assay <- CreateChromatinAssay(counts = data,sep = c("_", "_"), genome = genome, fragments = fragments, min.cells = 1, verbose = F)
metadata <- read.csv(file = metadata, header = TRUE, row.names = 1, sep="\t")
metadata <- metadata[ metadata$cell_called == "True",]
seurat <- CreateSeuratObject(counts = assay, assay = 'peaks', project = 'ATAC', meta.data = metadata)
Annotation(seurat) <- annotations     # add the gene information to the object

seurat <- NucleosomeSignal(object = seurat, verbose = F)

seurat$nucleosome_signal[ is.na(seurat$nucleosome_signal) ] <- 0    # NA to 0
seurat$nucleosome_signal[ !is.finite(seurat$nucleosome_signal) ] <- 0   # inf to 0

# head(Annotation(seurat))
# head(Fragments(seurat)[[1]])
seurat <- TSSEnrichment(seurat, fast = FALSE, verbose = F)

# QC plot ----
# nCount_peaks_TSS_enrichment
options(repr.plot.height = 5, repr.plot.width = 7)
DensityScatter(seurat, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(paste0(outdir, "/", name, "_nCount_peaks_TSS_enrichment.png"), width = 7, height = 5)

options(repr.plot.height = 5, repr.plot.width = 8)
seurat$high.tss <- ifelse(seurat$TSS.enrichment > 3, 'TSS_high', 'TSS_low')
TSSPlot(seurat, group.by = 'high.tss') + NoLegend()
ggsave(paste0(outdir, "/", name, "_high_TSS_enrichment.png"), width = 8, height = 5)

seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > 4, 'NS_high', 'NS_low')
table(seurat$nucleosome_group)

options(repr.plot.height = 5, repr.plot.width = 12)
VlnPlot(object = seurat, features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks', 'frac_peak'), pt.size = 0, ncol = 5)
ggsave(paste0(outdir, "/", name, "_QC.png"), width = 12, height = 5)

saveRDS(seurat, paste0(outdir, "/", name, ".rds"))
print("Done.")