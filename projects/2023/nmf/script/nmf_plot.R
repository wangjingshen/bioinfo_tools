options(warn = -1)
suppressMessages({
    library(argparser)
    library(Seurat)
    library(tidyverse)
    library(NMF)
})
options(warn = 1)

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--nmf_rds", help = "nmf rds")
argv <- add_argument(argv, "--outdir", help = "outdir to save fig. Default: . ")
argv <- add_argument(argv, "--height_consensusmap", help = "consensusmap height. Default: 10 ")
argv <- add_argument(argv, "--width_consensusmap", help = "consensusmap width. Default: 14 ")
argv <- add_argument(argv, "--height_heatmap", help = "heatmap height. Default: 5 ")
argv <- add_argument(argv, "--width_heatmap", help = "heatmap width. Default: 7 ")
argv <- add_argument(argv, "--topn", help = "topn features. Default: 5 ")

argv <- parse_args(argv)
argv$outdir <- ifelse(is.na(argv$outdir), ".", argv$outdir)
argv$height_consensusmap <- ifelse(is.na(argv$height_consensusmap), 20, as.numeric(argv$height_consensusmap))
argv$width_consensusmap <- ifelse(is.na(argv$width_consensusmap), 24, as.numeric(argv$width_consensusmap))
argv$height_heatmap <- ifelse(is.na(argv$height_heatmap), 6, as.numeric(argv$height_heatmap))
argv$width_heatmap <- ifelse(is.na(argv$width_heatmap), 7, as.numeric(argv$width_heatmap))
argv$topn <- ifelse(is.na(argv$topn), 5, as.numeric(argv$topn))

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive=T)
}

res <- readRDS(argv$nmf_rds)
#png( paste0(argv$outdir, "/consensusmap.png"), height = argv$height_consensusmap, width = argv$width_consensusmap)
#tiff( paste0(argv$outdir, "/consensusmap.tiff"))
#options(repr.plot.height = 10, repr.plot.width = 14)
tiff( paste0(argv$outdir, "/consensusmap.tiff"))
consensusmap(res, Colv=FALSE, Rowv=FALSE) 
dev.off()


print("nmf plot done.")