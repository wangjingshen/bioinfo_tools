options(warn = -1)
suppressMessages({
    library(argparser)
    library(Seurat)
    library(tidyverse)
    library(NMF)
})
options(warn = 1)

#
N_RUN <- 10
SEED <- 100

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--rank", help = "rank. Default: auto")
argv <- add_argument(argv, "--auto_max_rank", help = "auto max rank. Default: 10")
argv <- add_argument(argv, "--outdir", help = "outdir to save fig. Default: . ")
argv <- add_argument(argv, "--height_consensusmap", help = "consensusmap height. Default: 10 ")
argv <- add_argument(argv, "--width_consensusmap", help = "consensusmap width. Default: 14 ")
argv <- add_argument(argv, "--height_heatmap", help = "heatmap height. Default: 5 ")
argv <- add_argument(argv, "--width_heatmap", help = "heatmap width. Default: 7 ")
argv <- add_argument(argv, "--topn", help = "topn features. Default: 5 ")


argv <- parse_args(argv)

argv$rank <- ifelse(is.na(argv$rank), "auto", as.numeric(argv$rank))
argv$auto_max_rank <- ifelse(is.na(argv$auto_max_rank), 10, as.numeric(argv$auto_max_rank))
argv$outdir <- ifelse(is.na(argv$outdir), ".", argv$outdir)
argv$height_consensusmap <- ifelse(is.na(argv$height_consensusmap), 10, as.numeric(argv$height_consensusmap))
argv$width_consensusmap <- ifelse(is.na(argv$width_consensusmap), 14, as.numeric(argv$width_consensusmap))
argv$height_heatmap <- ifelse(is.na(argv$height_heatmap), 6, as.numeric(argv$height_heatmap))
argv$width_heatmap <- ifelse(is.na(argv$width_heatmap), 7, as.numeric(argv$width_heatmap))
argv$topn <- ifelse(is.na(argv$topn), 5, as.numeric(argv$topn))

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive=T)
}


# read data --
data_seurat <- readRDS(argv$rds)
data_seurat = CreateSeuratObject(counts = data_seurat@assays$RNA@counts, meta.data = data_seurat@meta.data)

data_seurat <- NormalizeData(data_seurat, verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(do.center = F, verbose = F)

# auto --
if(argv$rank == "auto"){
    res_auto <- nmf(data_seurat@assays$RNA@scale.data, rank=2:10, nrun=N_RUN, seed = SEED)
    rank <- which.min(diff(res_auto$measures$cophenetic)) + 1
}else{
    rank <- argv$rank
}

res <- nmf(data_seurat@assays$RNA@scale.data, rank=rank, nrun=N_RUN, seed = SEED)
saveRDS(res, paste0(argv$outdir,"/nmf_res.rds"))
write.table(rank, paste0(argv$outdir, "/best_rank.txt"))

# consensusmap --
tiff( paste0(argv$outdir, "/consensusmap.tiff"))
consensusmap(res, Colv=FALSE, Rowv=FALSE) 
dev.off()

# nmf features --
nmf_features_list <- lapply(extractFeatures(res), function(x) rownames(res)[x])
names(nmf_features_list) <- paste0("nmf_",1:length(nmf_features_list))
lapply(names(nmf_features_list), function(x){
    write.table(nmf_features_list[[x]], paste0(argv$outdir, "/", x, "_features.txt"), sep="\t", quote=F, row.names=F, col.names=F)
})

# nmf plot
nmf_features_topn <-  do.call(c,lapply(extractFeatures(res, argv$topn), function(x) rownames(res)[x]))
data_seurat_plot <- ScaleData(data_seurat, features = row.names(data_seurat@assays$RNA@counts), verbose = F)
DoHeatmap(object = data_seurat_plot, features = nmf_features_topn)
ggsave(paste0(argv$outdir, "/nmf_features_top", topn, "_heatmap.png"), height = argv$height_heatmap, width = argv$width_heatmap)


print("nmf done.")