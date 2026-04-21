suppressMessages({
    library(argparser)
    library(Seurat)
    library(tidyverse)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
})

# args ----
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--sweet_tsne_tag", help = "path of tsne_tag of sweet")
argv <- add_argument(argv, "--species", help = "species, human or mouse. Default: human")
#argv <- add_argument(argv, "--prefix", help = "(sample) prefix of barcode in rds")  # use sample_name
argv <- add_argument(argv, "--filter_cluster", help = "filter cluster, split by ','. Default: [Dd]oublet(s)")
argv <- add_argument(argv, "--max_cutoff", help = "max_cutoff for plot. Default: 10000")
argv <- add_argument(argv, "--analysis_cluster", help = "cluster for differential analysis. Default: The cluster with the highest number of cells")
argv <- add_argument(argv, "--percent", help = "percent of high(low) sweet tag. Default: 0.2")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: outdir")
argv <- add_argument(argv, "--width_go_high", help = "width of go res of high. Default: 10")
argv <- add_argument(argv, "--height_go_high", help = "height of go res of high. Default: 8")
argv <- add_argument(argv, "--width_go_low", help = "width of go res of low. Default: 10")
argv <- add_argument(argv, "--height_go_low", help = "height of go res of low. Default: 8")
argv <- parse_args(argv)

# arg default --
argv$species <- ifelse(is.na(argv$species), 'human', argv$species)
filter_cluster <- ifelse(is.na(argv$filter_cluster), c('Doublet','Doublets','doublet','doublets'), unlist(str_split(argv$filter_cluster, ',')))
argv$max_cutoff <- ifelse(is.na(argv$max_cutoff), 10000, as.numeric(argv$max_cutoff))
argv$percent <- ifelse(is.na(argv$percent), 0.2, as.numeric(argv$percent))
argv$outdir <- ifelse(is.na(argv$outdir), 'outdir', argv$outdir)
argv$width_go_high <- ifelse(is.na(argv$width_go_high), 10, as.numeric(argv$width_go_high))
argv$height_go_high <- ifelse(is.na(argv$height_go_high), 8, as.numeric(argv$height_go_high))
argv$width_go_low <- ifelse(is.na(argv$width_go_low), 10, as.numeric(argv$width_go_low))
argv$height_go_low <- ifelse(is.na(argv$height_go_low), 8, as.numeric(argv$height_go_low))

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)
}

# analysis ----
data_seurat <- readRDS(argv$rds)
sample_name <- unique(data_seurat$sample)

sweet_tag <- read.table(argv$sweet_tsne_tag, sep="\t", header = T, row.names = 1)
sweet_tag$barcode <- paste0(sample_name, "_", row.names(sweet_tag))
row.names(sweet_tag) <- sweet_tag$barcode

sweet_tag <- sweet_tag[colnames(data_seurat),]  # match 
if( length(intersect(sweet_tag$barcode, colnames(data_seurat)))!= ncol(data_seurat)){
    print("The Barcodes between rds and (subset) sweet are inconsistent, please check.")
    quit()
}else{
    data_seurat$tag_UMI <- sweet_tag$sweet_tag
    data_seurat <- subset(data_seurat, subset = cluster %in% setdiff(unique(data_seurat$cluster),filter_cluster))  # filter clusters

    # 1 sweet tag visualization --
    FeaturePlot(data_seurat, features ="tag_UMI", cols = c('lightgrey','#FF0000'), max.cutoff = argv$max_cutoff) +
        ggtitle(paste0(sample_name, ": tag_UMI")) +
        theme(plot.title = element_text(size=12))
    ggsave(paste0(argv$outdir, "/", sample_name, "_tag_UMI.pdf"), height = 5, width = 6)
    
    VlnPlot(data_seurat, features ='tag_UMI',group.by = "cluster", y.max = argv$max_cutoff, pt.size = 0) +
        NoLegend() +
        ggtitle(paste0(sample_name,": tag_UMI")) +
        theme(axis.title.x = element_blank(), plot.title = element_text(size=12))  
    ggsave(paste0(argv$outdir, "/", sample_name, "_tag_UMI_vlnplot.pdf"), height = 5, width = 6)

    function_set_sweet_group <- function(data, percent){
        group <- colnames(data)
        high_barcode <- names(data$tag_UMI[order(data$tag_UMI, decreasing = T)][1:floor(ncol(data)*percent)])
        group[ group %in% high_barcode] = "high"
        low_barcode <- names(data$tag_UMI[order(data$tag_UMI)][1:floor(ncol(data)*percent)])
        group[ group %in% low_barcode] = "low"
        group[ !group %in% c("high","low")] = "cm"
        return(group)
    }

    function_analysis <- function(cluster_name){
        data_analysis <- subset(data_seurat, subset = cluster == cluster_name)
        print(paste0("analyis cluster: ", cluster_name))
        print(paste0("cell number: ", ncol(data_analysis)))
        data_analysis$tag_UMI_group <- function_set_sweet_group(data_analysis, percent = argv$percent)
        data_analysis <- subset(data_analysis, subset = tag_UMI_group != "cm")  # filter cm
        data_analysis$tag_UMI_group <- factor(data_analysis$tag_UMI_group, levels = c('high','low'))
        print(paste0("cell number of high sweet: ", sum(data_analysis$tag_UMI_group == 'high')))
        print(paste0("cell number of low sweet: ", sum(data_analysis$tag_UMI_group == 'low')))

        # FindMarkers ----
        Idents(data_analysis) <- data_analysis$tag_UMI_group
        markers_df <- FindAllMarkers(data_analysis, only.pos = T, verbose = F)
        markers_df <- markers_df[markers_df$p_val_adj < 0.05,]
        markers_high = markers_df$gene[ markers_df$cluster=='high' & markers_df$p_val_adj < 0.05 ]
        markers_low = markers_df$gene[ markers_df$cluster=='low' & markers_df$p_val_adj < 0.05 ]
        write.table(markers_df[,c(7,1:6)], paste0(argv$outdir, "/", sample_name, "_", cluster_name, "_degs_table.tsv"), quote = F, sep="\t", row.names = F)

        # heatmap ----
        data_plot <- ScaleData(data_analysis, features = row.names(data_analysis@assays$RNA@data), verbose = F)    
        DoHeatmap(object = data_plot, 
                  features = c(markers_df %>% filter(cluster == "high") %>% top_n(n = 10, wt = avg_log2FC) %>% dplyr::select(gene) %>% unlist() %>% as.character(), 
                               markers_df %>% filter(cluster == "low") %>% top_n(n = 10, wt = avg_log2FC) %>% dplyr::select(gene) %>% unlist() %>% as.character()))
        ggsave(paste0(argv$outdir, "/", sample_name, "_", cluster_name, "_degs_heatmap.pdf"), width = 6, height = 5)

        # GO ----
        if(argv$species == 'human'){
            go_high <- enrichGO(gene = markers_high, OrgDb = "org.Hs.eg.db", keyType  = 'SYMBOL', ont="all")
            go_low <- enrichGO(gene = markers_low, OrgDb = "org.Hs.eg.db", keyType  = 'SYMBOL', ont="all")
        }else{
            go_high <- enrichGO(gene = markers_high, OrgDb = "org.Mm.eg.db", keyType  = 'SYMBOL', ont="all")
            go_low <- enrichGO(gene = markers_low, OrgDb = "org.Mm.eg.db", keyType  = 'SYMBOL', ont="all")
        }
        
        p1 <- dotplot(go_high, split="ONTOLOGY") +
            facet_grid(ONTOLOGY~., scale="free") +
            ggtitle(paste0(cluster_name, " :high")) +
            scale_y_discrete(labels = function(x) {str_wrap(x, width = 90)}) +
            scale_size()
        ggsave(paste0(argv$outdir, "/", sample_name, "_", cluster_name, "_go_high.pdf"), plot = p1, width = argv$width_go_high, height = argv$height_go_high)
        p2 <- dotplot(go_low, split="ONTOLOGY") +
            facet_grid(ONTOLOGY~., scale="free") + 
            ggtitle(paste0(cluster_name, ": low")) +
            scale_y_discrete(labels = function(x) {str_wrap(x, width = 90)}) + scale_size()
        ggsave(paste0(argv$outdir, "/", sample_name, "_", cluster_name, "_go_low.pdf"), plot = p2, width = argv$width_go_low, height = argv$height_go_low)
        saveRDS(list(go_high = go_high, go_low = go_low), paste0(argv$outdir, "/", sample_name, "_", cluster_name, "_go_res.rds"))  # save rds for graph adjustment
    }

    if(is.na(argv$analysis_cluster)){
        argv$analysis_cluster <- names(sort(table(data_seurat$cluster), decreasing = T))[1]
    }else{
        argv$analysis_cluster <- unlist(str_split(argv$analysis_cluster, ','))
    }
    lapply(argv$analysis_cluster, function_analysis)
    print('Done.')
}
