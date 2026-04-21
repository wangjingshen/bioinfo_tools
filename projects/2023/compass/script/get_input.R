options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--mtx", help = "matrix 10X")  # not recommended
argv <- add_argument(argv, "--subset", help = "subset")
argv <- add_argument(argv, "--subset_var", help = "subset variable")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir) 
}
argv$name <- ifelse(is.na(argv$name), "", paste0("_", argv$name))

if(!is.na(argv$rds)){
    data_seurat <- readRDS(argv$rds)
    data_counts <- data_seurat@assays$RNA@counts
    if(!is.na(argv$subset)){
        subset_cells <- data_seurat@meta.data[[argv$subset_var]] %in% unlist(strsplit(argv$subset, split = ","))
        data_counts <- data_counts[, subset_cells]
    }
    cpm <- apply(data_counts, 2, function(x){ x/sum(x)*10000 })
    write.table(cpm, paste0(argv$outdir, "/compass_input", argv$name, ".tsv"), quote = F, sep="\t")
}

if(!is.na(argv$mtx)){
    data_mtx <- Read10X(argv$mtx)
    cpm <- apply(data_mtx, 2, function(x){ x/sum(x)*10000 })
    write.table(cpm, paste0(argv$outdir, "/compass_input", argv$name, ".tsv"), quote = F, sep="\t")
}

print("Compass input done.")