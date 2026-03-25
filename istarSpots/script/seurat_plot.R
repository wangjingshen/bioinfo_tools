suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
# tab10
plot_col <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")

argv <- arg_parser('')
argv <- add_argument(argv, "--space_input", help = "path of space_input")
argv <- add_argument(argv, "--istar2spots", help = "istar2spots")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- parse_args(argv)

space_input <- argv$space_input
istar2spots <- argv$istar2spots
outdir <- argv$outdir

#
data_space <- Load10X_Spatial(space_input)

df <-read.table(istar2spots, sep=",", header = T)
df <- df[ match(colnames(data_space), df$barcode),]

if(identical(colnames(data_space), df$barcode)){
    data_space$istar_cluster = df$istar_cluster
    SpatialDimPlot(data_space, group.by = "istar_cluster", cols = plot_col, pt.size.factor = 1)
    ggsave(str_glue("{outdir}/istar_clusters.png"))
}else{
    print("barcode mismatch.")
    print(colnames(data_space)[1:5])
    print(df$barcode[1:5])
}

