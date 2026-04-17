# https://stuartlab.org/signac/articles/merging
# https://stuartlab.org/signac/articles/pbmc_multiomic

suppressWarnings(suppressMessages({
    library(Signac)
    library(Seurat)
    library(harmony)
    library(GenomicRanges)
    library(future)
    library(EnsDb.Mmusculus.v79)
    library(SingleR)
    library(tidyverse)
    library(argparser)
}))

plan("multicore", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2) # for 10 Gb RAM
set.seed(1)
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--fragments", help="fragments")
argv <- add_argument(argv,"--peak_bed", help="peak_bed")
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--mod", help="mod, 10X or self")
argv <- add_argument(argv,"--rm_batch", help="rm batch, yes or no. Default: no")
argv <- add_argument(argv,"--scRNA", help="scRNA rds")
argv <- add_argument(argv,"--species", help="species, human or mouse, Default: mouse")
argv <- add_argument(argv,"--metadata", help="metadata")
argv <- add_argument(argv,"--filter_bc", help="filter_bc")
argv <- add_argument(argv,"--resolution", help="resolution, Default: 0.8")
argv <- add_argument(argv,"--plot_region", help="plot_region, Default: chr14-99700000-99760000")
argv <- add_argument(argv,"--plot_cluster", help="plot_cluster, such as: Neurons")
argv <- add_argument(argv,"--plot_cluster_region", help="plot_cluster_region, such as: chr11-98279644-98379645")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- parse_args(argv)

fragments <- unlist(strsplit(argv$fragments, split = ","))
peak_bed <- unlist(strsplit(argv$peak_bed, split = ","))
sample <- unlist(strsplit(argv$sample, split = ","))
mod <- unlist(strsplit(argv$mod, split = ","))
rm_batch <- ifelse(is.na(argv$rm_batch), "no", argv$rm_batch)
scRNA <- argv$scRNA
species <- ifelse(is.na(argv$species), "mouse", argv$species)
metadata <- unlist(strsplit(argv$metadata, split = ","))
filter_bc <- unlist(strsplit(argv$filter_bc, split = ","))
resolution <- ifelse(is.na(argv$resolution), 0.8, as.numeric(argv$resolution))
plot_region <- ifelse(is.na(argv$plot_region), "chr14-99700000-99760000", argv$plot_region)
outdir <- argv$outdir

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}


# peak to granges
function_peak_to_granges <- function(peak_bed){
    # read in peak sets
    peak_bed <- read.table(peak_bed, col.names = c("chr", "start", "end"))
    # convert to genomic ranges
    granges_1 <- makeGRangesFromDataFrame(peak_bed)
    return(granges_1)
}

granges_list <- lapply(peak_bed, function_peak_to_granges)
names(granges_list) <- sample

# Create a unified set of peaks to quantify in each dataset
#combined.peaks <- reduce(x = c(granges_list$pipe_10X, granges_list$self_pipe))
combined.peaks <- granges_list[[1]]
if(length(granges_list)>1){
    for (i in 2:length(granges_list)){
        combined.peaks <- c(combined.peaks, granges_list[[i]])
    }
}
combined.peaks <- GenomicRanges::reduce(x = combined.peaks)

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# metadata
function_filter_bc_metadata <- function(x){
    if(x$mod == "10X"){
        metadata <- read.table(x$metadata, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ] # remove the first row
        filter_bc <- read.table(x$filter_bc)
        metadata_filter <- metadata[ filter_bc$V1,]
        return(metadata_filter)
    }
    if(x$mod == "self"){
        metadata <- read.table(x$metadata, stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = 1)
        metadata_filter <- metadata[ metadata$cell_called == "True",]
        return(metadata_filter)
    }
}

df_filter_bc_metadata <- data.frame(metadata = metadata,
                                    filter_bc = filter_bc,
                                    mod = mod)

metadata_filter_list <- lapply(split(df_filter_bc_metadata, seq(nrow(df_filter_bc_metadata))), function_filter_bc_metadata)
names(metadata_filter_list) <- sample

function_create_seurat <- function(x){
    # create fragment objects
    metadata <- metadata_filter_list[[x$sample]]
    frags <- CreateFragmentObject(x$fragments, cells = rownames(metadata))
    counts_1 <- FeatureMatrix( fragments = frags, features = combined.peaks, cells = rownames(metadata))
    
    assay_1 <- CreateChromatinAssay(counts_1, fragments = frags)
    seurat <- CreateSeuratObject(assay_1, assay = "ATAC", meta.data = metadata)
    seurat$sample <- x$sample
    return(seurat)
    #print(paste0(sample," done."))
}

df_create_seurat <- data.frame(fragments = fragments,
                               sample = sample)

#print(split(df_create_seurat, seq(nrow(df_create_seurat)))[[1]])
#test1 <- function_create_seurat(split(df_create_seurat, seq(nrow(df_create_seurat)))[[1]])
#print("1")
#test2 <- function_create_seurat(split(df_create_seurat, seq(nrow(df_create_seurat)))[[2]])
#print("2")

seurat_list <- lapply(split(df_create_seurat, seq(nrow(df_create_seurat))), function_create_seurat)

# merge all datasets, adding a cell ID to make sure cell names are unique
if(length(seurat_list)>1){
    combined <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample)
}else{
    combined <- seurat_list[[1]]
}
#combined[["ATAC"]]

if(rm_batch == "no"){
    combined <- RunTFIDF(combined, verbose = F)
    combined <- FindTopFeatures(combined, min.cutoff = 20, verbose = F)
    combined <- RunSVD(combined, verbose = F)
    combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi', verbose = F)
    combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30, verbose = F)
    combined <- FindClusters(object = combined, algorithm = 3, resolution = resolution, verbose = FALSE)
}else{
    combined <- RunTFIDF(combined)
    combined <- FindTopFeatures(combined, min.cutoff = 20)
    combined <- RunSVD(combined)
    combined <- NormalizeData(object = combined)
    combined <- ScaleData(object = combined)
    combined <- FindVariableFeatures(object = combined)
    genes.use<- head(HVFInfo(object = combined), 2000)
    combined <- RunPCA(object=combined, features = VariableFeatures(object = combined))
    combined <- RunHarmony(combined, group.by="sample" , plot_convergence = TRUE)
    combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
    combined <- FindClusters(combined, resolution = resolution, algorithm = 1)
    combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)
    #combined <- RunTSNE(object=combined, reduction = "harmony", dims.use=1:20, do.fast=TRUE, check_duplicates = FALSE)
}

## reference: https://stuartlab.org/signac/articles/pbmc_multiomic 
# extract gene annotations from EnsDb
if(species == "mouse"){
    suppressWarnings(annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79))
    # change to UCSC style since the data was mapped to hg19
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "mm10"
}
if(species == "human"){
    suppressWarnings(annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75))
    # change to UCSC style since the data was mapped to hg19
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "hg19"
}


# add the gene information to the object
Annotation(combined) <- annotations

# compute gene activities
gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(object = combined, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(combined$nCount_RNA))

# reference
reference <- readRDS(scRNA)
transfer.anchors <- FindTransferAnchors(reference = reference, query = combined, query.assay = "RNA", reduction = 'cca', verbose = F)

predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = reference$cluster, weight.reduction = combined[['lsi']], dims = 2:30, verbose = F)
combined <- AddMetaData(object = combined, metadata = predicted.labels)

options(repr.plot.width=7, repr.plot.height=5)
DimPlot(object = combined, group.by = "predicted.id", label = TRUE, repel = T)


options(repr.plot.height = 5, repr.plot.width = 7)
p <- DimPlot(object = combined, label = TRUE)
ggsave(paste0(outdir,"/cluster.pdf"), plot = p, height = 5, width = 7)
ggsave(paste0(outdir,"/cluster.png"), plot = p, height = 5, width = 7)

p <- DimPlot(combined, group.by = 'sample', pt.size = 0.1)
ggsave(paste0(outdir, "/sample.pdf"), plot = p, height = 5, width = 7)
ggsave(paste0(outdir, "/sample.png"), plot = p, height = 5, width = 7)

options(repr.plot.width = 8, repr.plot.height = 6)
p <- CoveragePlot(object = combined, group.by = 'sample', region = plot_region)
ggsave(paste0(outdir, "/CoveragePlot_", plot_region, ".pdf"), plot = p, height = 2 + 2*length(sample), width = 8)    #
ggsave(paste0(outdir, "/CoveragePlot_", plot_region, ".png"), plot = p, height = 2 + 2*length(sample), width = 8)    #

options(repr.plot.height = 6, repr.plot.width = 15)
plot1 <- DimPlot(object = reference, group.by = 'cluster', label = TRUE, repel = TRUE)  + ggtitle('scRNA-seq')
plot2 <- DimPlot(object = combined, group.by = 'predicted.id', label = TRUE, repel = TRUE) + ggtitle('scATAC-seq')
p <- plot1 + plot2
ggsave(paste0(outdir, "/anno.pdf"), plot = p, height = 6, width = 15)
ggsave(paste0(outdir, "/anno.png"), plot = p, height = 6, width = 15)

options(repr.plot.width=12, repr.plot.height = 5)
p <- DimPlot(combined, reduction = "umap", group.by = "predicted.id", label = T, repel = T, split.by = "sample") + NoLegend() + ggtitle("scATAC-seq")
ggsave(paste0(outdir, "/anno_split.pdf"), plot = p, height = 5, width = 3 + 3*length(sample))            #
ggsave(paste0(outdir, "/anno_split.png"), plot = p, height = 5, width = 3 + 3*length(sample))            #


if(!is.na(argv$plot_cluster)){
    plot_cluster <- argv$plot_cluster
    plot_cluster_region <- ifelse(is.na(argv$plot_cluster_region), "chr11-98279644-98379645", argv$plot_cluster_region)

    data_cluster <- subset(combined, subset = predicted.id == plot_cluster)
    options(repr.plot.height = 6, repr.plot.width = 9)
    p <- CoveragePlot(object = data_cluster, region = "chr11-98279644-98379645", group.by = 'sample', extend.upstream = 0, extend.downstream = 0)
    ggsave(paste0(outdir, "/CoveragePlot_", plot_cluster, "_", plot_cluster_region, ".pdf"), plot = p, height = 6, width = 9)
    ggsave(paste0(outdir, "/CoveragePlot_", plot_cluster, "_", plot_cluster_region, ".png"), plot = p, height = 6, width = 9)
    #CoveragePlot(object = data_neurons, region = "chr6-55631262-55731263", group.by = 'sample', extend.upstream = 0, extend.downstream = 0)
    #CoveragePlot(object = data_neurons, region = "chr4-45776826-45876827", group.by = 'sample', extend.upstream = 0, extend.downstream = 0)
}

saveRDS(combined, paste0(outdir, "/combined.rds"))
print("Done.")