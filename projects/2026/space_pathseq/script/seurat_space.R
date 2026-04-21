suppressWarnings(suppressMessages({
    library(Seurat)
    library(SeuratData)
    library(tidyverse)
    library(patchwork)
    library(dplyr)
    library(logger)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--pathseq_df", help = "pathseq df")
argv <- add_argument(argv, "--topn", default = 10, help="topn genus, default:10")
argv <- add_argument(argv, "--outdir", default = "outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

outdir <- argv$outdir

# space
data_seurat <- Load10X_Spatial("seurat_input/") %>%
    SCTransform(assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = F) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = F) 

# pathseq
pathseq_df <- read.table(argv$pathseq_df, sep="\t", header = T, row.names = 1)
data_seurat@meta.data <- cbind.data.frame(data_seurat@meta.data, t(pathseq_df))

stat_df <- data.frame(Genus = row.names(pathseq_df),
                      Number_of_UMIs = rowSums(pathseq_df),
                      Number_of_positive_UMIs = rowSums(pathseq_df > 0))
stat_df <- stat_df[ order(stat_df$Number_of_UMIs, decreasing = T),]

write.table(stat_df, str_glue("{outdir}/genus.tsv"), sep="\t", row.names = F, quote = F)

topn <- min(nrow(stat_df), argv$topn)

lapply(1:topn, function(x){
    plot_genus <- stat_df$Genus[x]
    SpatialFeaturePlot(object = data_seurat, features = plot_genus, alpha = c(0.1, 1))
    ggsave(str_glue("{outdir}/top{x}_{plot_genus}.png"), height = 6, width = 8)
})
