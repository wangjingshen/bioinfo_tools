suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(patchwork)
    library(argparser)
}))


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--spname", help = "spname")
argv <- add_argument(argv, "--biotype_gene_stats", help = "biotype_gene_stats")
argv <- add_argument(argv, "--ncRNA_type", default = "lncRNA,snoRNA,snRNA,miRNA,misc_RNA,scaRNA", help = "6 ncRNA type")
argv <- add_argument(argv, "--image_alpha", help = "image_alpha")
argv <- add_argument(argv, "--outdir", help="outdir")
argv <- parse_args(argv)

rds <- argv$rds
spname <- argv$spname
biotype_gene_stats <- argv$biotype_gene_stats
ncRNA <- unlist(strsplit(argv$ncRNA, split = ","))
image_alpha <- argv$image_alpha
outdir <- argv$outdir


function_test <- function(rds, biotype, subset_biotype){
    umi_matrix <- rds@assays$Spatial@counts
    subset_genes <- intersect(row.names(umi_matrix), biotype$gene_symbol[ biotype$gene_biotype == subset_biotype])
    rds@meta.data[, paste0(subset_biotype, "_total_UMI")] <- colSums(umi_matrix[subset_genes, ])
    rds@meta.data[, paste0(subset_biotype, "_log2UMI")] <- log2(rds@meta.data[, paste0(subset_biotype, "_total_UMI")] + 1)
    return(rds)
}

## read
data <- readRDS(rds)
biotype <- read.table(biotype_gene_stats, sep=",", header = T)

for (i in ncRNA){
    data <- function_test(data, biotype, i)
}

ncRNA_info <- list(
    list(feat = "lncRNA_total_UMI", low = "white", high = "#1f78b4"),
    list(feat = "snoRNA_total_UMI", low = "white", high = "#238b45"),
    list(feat = "snRNA_total_UMI",  low = "white", high = "#ffb300"),
    list(feat = "miRNA_total_UMI", low = "white", high = "#7a0177"),
    list(feat = "misc_RNA_total_UMI", low = "white", high = "#8c510a"),
    list(feat = "scaRNA_total_UMI", low = "white", high = "#e6550d")
)
plot_list <- list()
for (item in ncRNA_info) {
    p <- SpatialFeaturePlot(data, features = item$feat, alpha = c(0.1, 1),pt.size.factor = 1.6,image.alpha = image_alpha) +
            ggtitle(item$feat)+
            scale_fill_gradientn(colors = c(item$low, item$high)) + #自定义渐变色
            theme(legend.position = "right", plot.title=element_text(hjust = 0.5), legend.title=element_blank())
    plot_list[[item$feat]] <- p
}
final_plot <- wrap_plots(plot_list, ncol = 3)
ggsave(str_glue("plot/featureplot_ncRNA_total_UMI_{spname}.png"), width = 15, height = 9)
ggsave(str_glue("plot/featureplot_ncRNA_total_UMI_{spname}.pdf"), width = 15, height = 9)


ncRNA_info <- list(
    list(feat = "lncRNA_log2UMI", low = "white", high = "#1f78b4"),
    list(feat = "snoRNA_log2UMI", low = "white", high = "#238b45"),
    list(feat = "snRNA_log2UMI",  low = "white", high = "#ffb300"),
    list(feat = "miRNA_log2UMI", low = "white", high = "#7a0177"),
    list(feat = "misc_RNA_log2UMI", low = "white", high = "#8c510a"),
    list(feat = "scaRNA_log2UMI", low = "white", high = "#e6550d")
)
plot_list <- list()
for (item in ncRNA_info) {
    p <- SpatialFeaturePlot(data, features = item$feat, alpha = c(0.1, 1),pt.size.factor = 1.6,image.alpha = image_alpha) +
            ggtitle(item$feat)+
            scale_fill_gradientn(colors = c(item$low, item$high)) + #自定义渐变色
            theme(legend.position = "right", plot.title=element_text(hjust = 0.5), legend.title=element_blank())
    plot_list[[item$feat]] <- p
}
final_plot <- wrap_plots(plot_list, ncol = 3)
ggsave(str_glue("plot/featureplot_ncRNA_log2UMI_{spname}.png"), width = 15, height = 9)
ggsave(str_glue("plot/featureplot_ncRNA_log2UMI_{spname}.pdf"), width = 15, height = 9)
