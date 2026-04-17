options(warn = -1)    # off warnings 
suppressMessages({
    library(argparser)
    library(tidyverse)})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: HBV_stat")
argv <- parse_args(argv)

# 
data_seurat <- readRDS(argv$rds)
write.table(data_seurat@assays$RNA@data, paste0(argv$outdir,"mtx.txt"), sep="\t", quote = F, row.names = F)

print('Get mtx from rds done.')