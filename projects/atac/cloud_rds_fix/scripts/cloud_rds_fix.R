suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

argv <- arg_parser('')
argv <- add_argument(argv,"--rds",help="Seurat rds")
argv <- add_argument(argv,"--set_group_by_sample", help="set_group_by_sample, default:F")
argv <- add_argument(argv,"--outdir", help="the output dir. default: outdir")
argv <- parse_args(argv)

#read args
rds <- argv$rds
set_group_by_sample <- ifelse(is.na(argv$set_group_by_sample), "F", argv$set_group_by_sample)
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

data_seurat <- readRDS(rds)
print(data_seurat@meta.data[1:2,])
# fix
data_seurat$sample <- data_seurat$`Sample ID`
if(set_group_by_sample == "T"){
    data_seurat$group <- data_seurat$sample
}
data_seurat$cluster <- data_seurat$annot_full
Idents(data_seurat) <- data_seurat$cluster

saveRDS(data_seurat, str_glue("{outdir}/data.rds"))