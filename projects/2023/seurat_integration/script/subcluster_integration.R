options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
})
options(warn = 1)

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--batch_var", help = "batch variable name")
argv <- add_argument(argv, "--celltype_var", help = "cell type variable name")
argv <- add_argument(argv, "--subset_celltype", help = "subset cell type for subclustering")
argv <- add_argument(argv, "--nfeatures", help = "nfeatures, default: 2000")
argv <- add_argument(argv, "--npcs", help = "npcs, default: 20")
argv <- add_argument(argv, "--resolution", help = "resolution, default: 0.6")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

nfeatures <- ifelse(is.na(argv$nfeatures), 2000, argv$nfeatures)
npcs <- ifelse(is.na(argv$npcs), 20, argv$npcs)
resolution <- ifelse(is.na(argv$resolution), 0.6, argv$resolution)
outdir <- ifelse(is.na(argv$outdir), './', argv$outdir)

# read data --
data_seurat <- readRDS(argv$rds)
subset_celltype <- unlist(strsplit(argv$subset_celltype, split = ","))
data_seurat <- subset(data_seurat, subset = celltype %in% subset_celltype)   # subset for subclustering

data_list <- SplitObject(data_seurat, split.by = argv$batch_var)
data_seurat_list <- lapply(1:length(data_list), function(i){
    data <- CreateSeuratObject(data_list[[i]]@assays$RNA@counts)
    data$batch <- data_list[[i]][[argv$batch_var]]
    if(!is.na(argv$celltype_var)){
        data$celltype <- data_list[[i]][[argv$celltype_var]]
    }
    data <- NormalizeData(data)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures)
    return(data)
})

# integration --
anchor_features <- SelectIntegrationFeatures(object.list = data_seurat_list)
anchors <- FindIntegrationAnchors(object.list = data_seurat_list, anchor.features = anchor_features)
data_integrate <- IntegrateData(anchorset = anchors)
DefaultAssay(data_integrate) <- "integrated"

# clustering --
data_integrate <- data_integrate %>% 
  ScaleData() %>%
  RunPCA(verbose = F) %>%
  FindNeighbors(dims = 1:npcs, reduction = "pca", verbose = F) %>%
  FindClusters(resolution = resolution, verbose = F) %>%
  RunTSNE(dims=1:npcs, do.fast = TRUE, check_duplicates = FALSE, verbose = F) %>%
  RunUMAP(dims=1:npcs, reduction = "pca", verbose = F)

if(is.na(argv$name)){
    saveRDS(data_integrate, paste0(outdir, "/data_subclustering_integrate.rds"))
}else{
    saveRDS(data_integrate, paste0(outdir, "/", argv$name, ".rds"))
}