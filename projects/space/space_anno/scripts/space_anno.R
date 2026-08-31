# https://satijalab.org/seurat/articles/spatial_vignette

suppressWarnings(suppressMessages({
    library(Seurat)
    library(SeuratData)
    library(tidyverse)
    library(patchwork)
    library(dplyr)
    library(argparser)
}))

color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8",
"#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F",
"#000101","OrangeRed","SlateBlue","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink","Red","#4682B4",
"#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00",
"#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE",
"#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD",
"#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF",
"#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4",
"#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F",
"#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347",
"#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF",
"#1E90FF","#191970","#E8E8E8","#FFDAB9")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--space_input", help = "space input")
argv <- add_argument(argv, "--sc", help = "sc anno rds")
argv <- add_argument(argv, "--pt_size", default = 1.6, help = "pt_size, Default:1.6")
argv <- add_argument(argv, "--score_filter", help = "score filter")
argv <- add_argument(argv, "--filter_cluster", help = "filter cluster in sc")
argv <- add_argument(argv, "--resolution", default = 0.5, help = "resolution")
argv <- add_argument(argv, "--name", help="name")
argv <- parse_args(argv)

space_input <- argv$space_input
sc <- argv$sc
pt_size <- as.numeric(argv$pt_size)
score_filter <- as.numeric(argv$score_filter)
filter_cluster <- unlist(str_split(argv$filter_cluster, ","))
resolution <- argv$resolution
name <- argv$name
if(!dir.exists(str_glue("{name}/plot"))){
    dir.create(str_glue("{name}/plot"), recursive = TRUE)
}

# read space
data_space <- Load10X_Spatial(space_input)
# fix 10X demo data
#print(data_space@images$slice1@coordinates[1:2,])
#print(class(data_space@images$slice1@coordinates$imagerow))  # character -> numeric
data_space@images$slice1@coordinates$imagerow <- as.numeric(data_space@images$slice1@coordinates$imagerow)
data_space@images$slice1@coordinates$imagecol <- as.numeric(data_space@images$slice1@coordinates$imagecol)

# filter 0 umi
data_space$total_umi <- colSums(data_space@assays$Spatial@counts)
data_space <- subset(data_space, cells = colnames(data_space)[data_space$total_umi > 0])

# space SCT
data_space <- data_space %>%
    SCTransform(assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = F) %>%
    FindClusters(resolution = resolution, verbose = FALSE)

# plot
SpatialDimPlot(data_space, group.by = "seurat_clusters", cols = color_protocol, pt.size.factor = pt_size)
ggsave(str_glue("{name}/plot/space_seurat_clusters.png"))
ggsave(str_glue("{name}/plot/space_seurat_clusters.pdf"))

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
ggsave(str_glue("{name}/plot/sc_clusters.png"), height = 6, width = 10)
ggsave(str_glue("{name}/plot/sc_clusters.pdf"), height = 6, width = 10)

anchors <- FindTransferAnchors(reference = data_sc, query = data_space, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = data_sc$cluster, prediction.assay = TRUE, 
                                  weight.reduction = data_space[["pca"]], dims = 1:30)

data_space[["predictions"]] <- predictions.assay

DefaultAssay(data_space) <- "predictions"
data_space$cluster <- GetTransferPredictions(data_space, assay = "predictions", slot = "data", score.filter = score_filter)

Idents(data_space) <- data_space$cluster  # "cluster"
cell_types_highlight <- CellsByIdentities(object = data_space, idents = unique(data_space$cluster))
SpatialDimPlot(data_space, cells.highlight = cell_types_highlight, facet.highlight = TRUE, ncol = ceiling(length(unique(data_space$cluster))/2))
per_width = 4
ggsave(str_glue("{name}/plot/space_cluster_split.png"), width = 4*length(unique(data_space$cluster)), height = 9)
ggsave(str_glue("{name}/plot/space_cluster_split.pdf"), width = 4*length(unique(data_space$cluster)), height = 9)

SpatialDimPlot(data_space, group.by = "cluster", cols = color_protocol, pt.size.factor = pt_size)
ggsave(str_glue("{name}/plot/space_cluster.png"), width = 8, height = 6)
ggsave(str_glue("{name}/plot/space_cluster.pdf"), width = 8, height = 6)

#
SpatialFeaturePlot(data_space, features = unique(data_space$cluster), alpha = c(0.1, 1), ncol = ceiling(length(unique(data_space$cluster))/2))
ggsave(str_glue("{name}/plot/space_predictions_score.png"), width = 4*length(unique(data_space$cluster)), height = 9)
ggsave(str_glue("{name}/plot/space_predictions_score.pdf"), width = 4*length(unique(data_space$cluster)), height = 9)


saveRDS(data_space, str_glue("{name}/data_space.rds"))