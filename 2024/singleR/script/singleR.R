
suppressWarnings(suppressMessages({
    library(argparser)
    library(Seurat)
    library(SingleR)
    library(tidyverse)
}))

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds")
argv <- add_argument(argv,"--species", help="species")
argv <- add_argument(argv,"--outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

rds <- argv$rds
species <- ifelse(is.na(argv$species), "human", argv$species)
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

function_singleR <- function(data, species){
    if(species == "human"){
        ref = readRDS("/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/HumanPrimaryCellAtlasData.rds")
    }
    if(species == "mouse"){
        ref = readRDS("/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/MouseRNAseqData.rds")
    }
    anno <- SingleR(test = data@assays$RNA@data, ref = ref, 
                    clusters = unlist(data$seurat_clusters),
                    assay.type.test=1, labels = ref$label.main)
    cell_type_singleR <- plyr::mapvalues(x = data$seurat_clusters, from = row.names(anno), to = anno$labels)
    return(cell_type_singleR)
}
data_seurat <- readRDS(argv$rds)
data_seurat$cluster_singleR <- function_singleR(data_seurat, species)
data_seurat$cluster_singleR <- factor(as.character(data_seurat$cluster_singleR), levels= sort(unique(as.character(data_seurat$cluster_singleR))))

# plot
p1 <- DimPlot(data_seurat,group.by = "cluster_singleR",label = T,repel = T)
ggsave(paste0(outdir,"/seurat_cluster_singleR.png"), plot=p1, height=6.5, width=8)
ggsave(paste0(outdir,"/seurat_cluster_singleR.pdf"), plot=p1, height=6.5, width=8)

# cellular_composition
p2 <- data_seurat@meta.data %>%
    ggplot(aes(x=sample,fill=cluster_singleR))+
        geom_bar(color="black",position="fill",width = 0.5)+
        theme_bw()+
        theme(legend.title = element_blank(),axis.title = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
ggsave(paste0(outdir, "/cellular_composition.png"), plot= p2, width = 4, height = 5)
ggsave(paste0(outdir, "/cellular_composition.pdf"), plot= p2, width = 4, height = 5)

# saveRDS
saveRDS(data_seurat, paste0(outdir, "/seurat_singleR.rds"))