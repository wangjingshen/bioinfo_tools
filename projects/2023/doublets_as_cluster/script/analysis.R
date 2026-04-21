suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "06 analysis rds")
argv <- add_argument(argv, "--mix_species", help = "mix species from doublets")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = T)
}

data_seurat <- readRDS(argv$rds)
doublets <- read.table(argv$mix_species, header = T)
doublets[1,1] <- colnames(data_seurat)[1]
if(identical(doublets$cell_id, colnames(data_seurat))){
    print("The barcode of doublets and barcode of rds are the same.")
    data_seurat$cluster <- doublets$Species
    saveRDS(data_seurat, paste0(argv$outdir, "/", argv$name,".rds"))
}else{
    print("The barcode of doublets and barcode of rds are different.")
}
