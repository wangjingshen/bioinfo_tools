suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--space_input", help = "space input")
argv <- add_argument(argv, "--tag_umi", help = "space input")
argv <- add_argument(argv, "--outdir", default = "outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

space_input <- argv$space_input
tag_umi <- argv$tag_umi
outdir <- argv$outdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}


data_space <- Load10X_Spatial(space_input)
tag_umi <- read.table(tag_umi, header=T)

tryCatch({
    if (!identical(colnames(data_space), row.names(tag_umi))) {
        stop("barcode_mismatch")
    }
}, error = function(e) {
    print("The barcodes do not match, please check.")
    quit()
})


data_space$tag <- tag_umi[,1]

data_space$log2_nCount_Spatial <- log2(data_space$nCount_Spatial + 1)
data_space$log2_tag <- log2(data_space$tag + 1)


data_space <- data_space %>%
    SCTransform(assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = F) %>%
    FindClusters(verbose = FALSE)

plot_vln <- function(feat) {
    p <- VlnPlot(data_space, features = feat, group.by = "seurat_clusters", pt.size = 0)
    ggsave(str_glue("{outdir}/VlnPlot_{feat}.png"), plot = p)
    p
}
feature_list <- c("nCount_Spatial", "tag", "log2_nCount_Spatial", "log2_tag")
p_list <- map(feature_list, plot_vln)