# https://satijalab.org/seurat/articles/spatial_vignette

suppressWarnings(suppressMessages({
    library(Seurat)
    library(SeuratData)
    library(tidyverse)
    library(patchwork)
    library(dplyr)
    library(argparser)
}))

color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--space_input", help = "space input")
argv <- add_argument(argv, "--sc", help = "sc anno rds")
argv <- add_argument(argv, "--score_filter", help = "score filter")
argv <- add_argument(argv, "--filter_cluster", help = "filter cluster in sc")
argv <- add_argument(argv, "--outdir", default = "outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

space_input <- argv$space_input
sc <- argv$sc
score_filter <- as.numeric(argv$score_filter)
filter_cluster <- unlist(str_split(argv$filter_cluster, ","))
outdir <- argv$outdir
if(!dir.exists(outdir)){
    dir.create(str_glue("{outdir}/plot"), recursive = TRUE)
}

# read space
data_space <- Load10X_Spatial(space_input)
# fix 10X demo data
#print(data_space@images$slice1@coordinates[1:2,])
#print(class(data_space@images$slice1@coordinates$imagerow))  # character
data_space@images$slice1@coordinates$imagerow <- as.numeric(data_space@images$slice1@coordinates$imagerow)
data_space@images$slice1@coordinates$imagecol <- as.numeric(data_space@images$slice1@coordinates$imagecol)
#print(class(data_space@images$slice1@coordinates$imagerow))  # numeric

# space SCT
data_space <- data_space %>%
    SCTransform(assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = F) %>%
    FindClusters(verbose = FALSE)

# plot
SpatialDimPlot(data_space, group.by = "seurat_clusters", cols = color_protocol)
ggsave(str_glue("{outdir}/plot/space_seurat_clusters.png"))

# read sc
data_sc <- readRDS(sc)
if("annot_full" %in% colnames(data_sc@meta.data)){
    data_sc$cluster <- gsub(" ", "_", data_sc$annot_full)  # sgr cloud annotation
}
data_sc <- subset(data_sc, subset = cluster %in% setdiff(unique(data_sc$cluster), filter_cluster))

## sc SCT
data_sc <- SCTransform(data_sc, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
DimPlot(data_sc, group.by = "cluster", cols = color_protocol, label = TRUE)
ggsave(str_glue("{outdir}/plot/sc_clusters.png"), height = 6, width = 10)

anchors <- FindTransferAnchors(reference = data_sc, query = data_space, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = data_sc$cluster, prediction.assay = TRUE, 
                                  weight.reduction = data_space[["pca"]], dims = 1:30)

data_space[["predictions"]] <- predictions.assay

DefaultAssay(data_space) <- "predictions"
data_space$predicted_id <- GetTransferPredictions(data_space, assay = "predictions", slot = "data", score.filter = score_filter)

Idents(data_space) <- data_space$predicted_id  # "predicted_id"
cell_types_highlight <- CellsByIdentities(object = data_space, idents = unique(data_space$predicted_id))
SpatialDimPlot(data_space, cells.highlight = cell_types_highlight, facet.highlight = TRUE, ncol = ceiling(length(unique(data_space$predicted_id))/2))
per_width = 4
ggsave(str_glue("{outdir}/plot/space_predicted_id_split.png"), width = 4*length(unique(data_space$predicted_id)), height = 9)

SpatialDimPlot(data_space, group.by = "predicted_id", cols = color_protocol)
ggsave(str_glue("{outdir}/plot/space_predicted_id.png"), width = 8, height = 6)

#
SpatialFeaturePlot(data_space, features = unique(data_space$predicted_id), alpha = c(0.1, 1), ncol = ceiling(length(unique(data_space$predicted_id))/2))
ggsave(str_glue("{outdir}/plot/space_predictions_score.png"), width = 4*length(unique(data_space$predicted_id)), height = 9)

saveRDS(data_space, str_glue("{outdir}/data_space.rds"))

# use max (test)
#predictions.res <- predictions.assay@data
#rownames(predictions.res) <- make.names(rownames(predictions.res))
#cell_types <- rownames(predictions.res)
#for(cell_type in cell_types){
#    data_space <- AddMetaData(object = data_space, metadata = predictions.res[cell_type,], col.name = cell_type)
#}
#data_space$clusters <- apply(predictions.res, 2, function(x){
#    rownames(predictions.res)[which.max(x)]
#})
#SpatialDimPlot(data_space, group.by = "clusters", cols = color_protocol)
#ggsave(str_glue("{outdir}/space_clusters.png"))

