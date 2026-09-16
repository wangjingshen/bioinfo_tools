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
argv <- add_argument(argv, "--space_input", help = "space input")
argv <- add_argument(argv, "--virus_df", help = "virus df")
argv <- add_argument(argv, "--name", help="name")
argv <- parse_args(argv)

space_input <- argv$space_input
virus_df <- argv$virus_df
name <- argv$name

# space
data_seurat <- Load10X_Spatial(space_input) 
# virus
virus_df <- read.table(argv$virus_df, sep=",", header = T, row.names = 1)
identical(colnames(data_seurat), row.names(virus_df))

data_seurat$virus_UMI <- virus_df$sum_UMI

print(median(data_seurat$virus_UMI))
print(mean(data_seurat$virus_UMI))

SpatialFeaturePlot(object = data_seurat, features = "virus_UMI", alpha = c(0.1, 1), image.alpha = 0.3) +
    scale_fill_gradient(low = "white", high = "#00A23F") +
    theme(legend.position = "right")
ggsave(str_glue("{name}/07.analysis_virus/{name}_virus_UMI.png"), height = 6, width = 8)
