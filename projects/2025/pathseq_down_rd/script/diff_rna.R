#

suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--genus_analysis", help = "genus to split group, total_genus or genus")
#argv <- add_argument(argv, "--umi_threshold", help = "umi threshold, total_genus:2, genus:0")
argv <- add_argument(argv, "--outdir", help = "output dir, default: outdir")
argv <- parse_args(argv)

rds <- argv$rds
genus_analysis <- ifelse(is.na(argv$genus_analysis), "total_genus", argv$genus_analysis)
umi_threshold <- ifelse(genus_analysis == "total_genus", 2, 0)
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(str_glue("{outdir}/02.diff_rna/{genus_analysis}"))){
    dir.create(str_glue("{outdir}/02.diff_rna/{genus_analysis}"), recursive = T)
}

data_seurat <- readRDS(rds)
data_seurat$genus_group <- factor(ifelse(data_seurat@meta.data[,genus_analysis] >0, "genus_positive", "genus_negative"), levels = c("genus_positive", "genus_negative"))
print(paste0("genus: ",genus_analysis))
print(table(data_seurat$cluster,data_seurat$genus_group))

function_diff <- function(subcluster){
    data_sub <- subset(data_seurat, subset = cluster == subcluster)
    if(length(unique(data_sub$genus_group))>1){    # 1: the cluster no genus detect
        Idents(data_sub) <- data_sub$genus_group
        markers <- FindMarkers(data_sub, ident.1 = "genus_positive", ident.2 = "genus_negative", verbose = F)
        markers$gene <- row.names(markers)
        colnames(markers)[2] <- "avg_logFC"
        print(subcluster)
        print(table(markers$avg_logFC>0))
        write.table(markers[,c(6,1:5)], str_glue("{outdir}/02.diff_rna/{genus_analysis}/degs_{subcluster}_{genus_analysis}_positive_vs_negative.tsv"), sep="\t", row.names = F, quote=F)
    }
    if(nrow(markers)>2){
        data_sub <- ScaleData(data_sub, features = row.names(data_sub@assays$RNA@counts), verbose = F)
        markers %>%
            top_n(n = 10, wt = avg_logFC) -> markers_top10
        markers %>%
            top_n(n = 10, wt = - avg_logFC) -> markers_bottom10
        DoHeatmap(data_sub, features = c(markers_top10$gene, markers_bottom10$gene), group.by = "genus_group", label = F)
        ggsave(str_glue("{outdir}/02.diff_rna/{genus_analysis}/heatmap_{subcluster}_{genus_analysis}.png"))
        ggsave(str_glue("{outdir}/02.diff_rna/{genus_analysis}/heatmap_{subcluster}_{genus_analysis}.pdf"))
    }
    return(markers)
}
markers <- lapply(levels(data_seurat$cluster), function_diff)
names(markers) <- levels(data_seurat$cluster)
cat("Diff done. \n")