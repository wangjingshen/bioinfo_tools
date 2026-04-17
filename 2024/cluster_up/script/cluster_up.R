suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)
}

data_seurat <- readRDS(argv$rds)
print(table(Idents(data_seurat)))
# diff
markers <- FindAllMarkers(data_seurat, only.pos = T,verbose = F)
markers_split <- split(markers, markers$cluster)
sapply(names(markers_split),function(x){
    write.table(markers_split[[x]][,c(7,1:6)], paste0(argv$outdir, "/", x, ".xls"), sep="\t", quote = F, row.names = F)
})

diffgenes_list <- data.frame(cluster_name = names(markers_split),
                             diffgenes_path = paste0(argv$outdir, "/", names(markers_split), ".xls"))
write.table(diffgenes_list, "diffgenes_list", sep="\t", row.names = F, col.names = F, quote = F)

print("Done.")