suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(enrichplot)
    library(patchwork)
    library(dplyr)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--spname", help = "spname")
argv <- add_argument(argv, "--idents", default ="seurat_clusters", help = "idents, default:seurat_clusters")
argv <- add_argument(argv, "--analysis_cluster", help="analysis_cluster")
argv <- add_argument(argv, "--compare_cluster", help="compare_cluster")
argv <- add_argument(argv, "--outdir", help="outdir")
argv <- parse_args(argv)

rds <- argv$rds
spname <- argv$spname
idents<- argv$idents
analysis_cluster <- argv$analysis_cluster
compare_cluster <- argv$compare_cluster
outdir <- argv$outdir


#
data <- readRDS(rds)
DefaultAssay(data) <- "SCT"

# dotplot
Idents(data) <- data@meta.data[,idents]
marker_res <- FindMarkers(data, ident.1 = analysis_cluster, ident.2 = compare_cluster, only.pos = FALSE, 
                          min.pct = 0.05, logfc.threshold = 0.25, test.use = "wilcox")
    
gene_all <- marker_res %>%
    filter(p_val_adj < 0.05) %>%
    group_by(direction = if_else(avg_log2FC > 0, "up", "down")) %>%
    arrange(direction, desc(avg_log2FC)) %>%
    slice_head(n = 10) %>%
    pull(rowname)

DotPlot(data, features = gene_all, idents = c(analysis_cluster, compare_cluster)) +
    RotatedAxis() +
    theme_test() +
    labs(x = "Genes", y = "Cell types") + 
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust =1))
ggsave(str_glue("{outdir}/dotplot_{analysis_cluster}vs{compare_cluster}_{spname}.pdf"), width = 12, height = 4, dpi = 300)
ggsave(str_glue("{outdir}/dotplot_{analysis_cluster}vs{compare_cluster}_{spname}.png"), width = 12, height = 4, dpi = 300)

# diff
set.seed(123)
data <- subset(data, downsample = 200)
all_markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25, test.use = "wilcox")

top_genes <- all_markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(n = n_top, order_by = avg_log2FC, with_ties = FALSE) %>%
    pull(gene)

p <- DoHeatmap(object = data, features = top10_genes, group.by = "seurat_clusters", assay = "SCT", size = 3) + NoLegend()
ggsave(str_glue("{outdir}/heatmap_markers_{spname}.png"), p, width = 10, height = 12)
ggsave(str_glue("{outdir}/heatmap_markers_{spname}.pdf"), p, width = 10, height = 12)
    
marker_gene_sym <- all_markers %>%
    filter(cluster == analysis_cluster, p_val_adj < 0.05) %>%
    pull(gene)

go_all_df <- tibble(ont = c("BP","CC","MF")) %>%
    mutate(go_res = purrr::map(ont, ~enrichGO(gene = marker_gene_sym, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                                              ont = .x, pAdjustMethod = "fdr", qvalueCutoff = 0.05))) %>%
    mutate(go_df = purrr::map(go_res, ~as.data.frame(.x))) %>%
    unnest(go_df) %>%
    select(-go_res)

go_top <- go_all_df %>%
    group_by(ont) %>%
    slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    separate(GeneRatio, c("num","den"), sep="/", convert = TRUE) %>%
    mutate(GeneRatio = num / den) %>%
    group_by(ont) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    mutate(Description = factor(Description, levels = unique(Description))) %>%
    ungroup()

p_go <- ggplot(go_top, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    facet_grid(ont ~ ., scale = "free_y") +
    scale_color_gradient(low = "#4444ff", high = "#ff2222") +
    scale_size_continuous(range = c(2, 8)) +
    scale_x_continuous(breaks = seq(0, max(go_top$GeneRatio, na.rm = TRUE), by = 0.05)) +
    labs(x="GeneRatio", y="Description", color="-log10(p.adjust)", size="Count") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 15),
          strip.text.y.right = element_text(size = 15, angle = 0),
          strip.placement = "right",
          panel.grid = element_line(color = "grey90"),
          legend.position = "right"
    )

ggsave(str_glue("{outdir}/GO_cluster{analysis_cluster}_{spname}.pdf"), plot = p_go, width=13, height=9)
ggsave(str_glue("{outdir}/GO_cluster{analysis_cluster}_{spname}.png"), plot = p_go, width=13, height=9)
