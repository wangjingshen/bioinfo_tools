### v2 for lims
# fj version 1.5.1; celescope1.5.1b0

suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(patchwork)})

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--version", help = "version. v1 or v2. v2 for lims. Default: v2")
argv <- parse_args(argv)

# default value -- 
if(is.na(argv$version)){
    argv$version ="v2"
}
# version = argv$version # version is R function, not set, use argv$version

# read rds -- 
data_seurat <- readRDS(argv$rds)
if(argv$version == "v1"){
    print(table(data_seurat$orig.ident))
}
if(argv$version == "v2"){
    print(table(data_seurat$sample))
}

# end --
cat("------------\n")
cat("get prefix done.\n")