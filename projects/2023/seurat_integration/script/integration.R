options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
})
options(warn = 1)

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--mode", help = "mode, matrix or rds")
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--batch_var", help = "batch variable, supoort mode: rds")
argv <- add_argument(argv, "--matrix_10X", help = "10X matrix, split by ,")
argv <- add_argument(argv, "--batch_name", help = "batch name, split by, support mode: matrix")
argv <- add_argument(argv, "--celltype_var", help = "celltype variable, for test pancreas data")
argv <- add_argument(argv, "--nfeatures", help = "nfeatures, default: 2000")
argv <- add_argument(argv, "--npcs", help = "npcs, default: 20")
argv <- add_argument(argv, "--resolution", help = "resolution, default: 0.6")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

nfeatures <- ifelse(is.na(argv$nfeatures), 2000, argv$nfeatures)
npcs <- ifelse(is.na(argv$npcs), 20, argv$npcs)
resolution <- ifelse(is.na(argv$resolution), 0.6, argv$resolution)
outdir <- ifelse(is.na(argv$outdir), './', argv$outdir)

#if(is.na(argv$npcs)){
#    npcs <- 20
#}else{
#    npcs <- argv$npcs
#}
#if(is.na(argv$resolution)){
#    resolution <- 0.6
#}else{
#    resolution <- argv$resolution
#}

if(argv$mode == "matrix"){
    matrix_10X <- unlist(strsplit(argv$matrix_10X, split = ","))
    batch_name <- unlist(strsplit(argv$batch_name, split = ","))
    data_count <- lapply(matrix_10X, Read10X)
    data_seurat_list <- lapply(1:length(data_count), function(i){
        data <- CreateSeuratObject(data_count[[i]])
        data$batch <- batch_name[[i]]
        data <- NormalizeData(data)
        data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        return(data)
    })
}

if(argv$mode == "rds"){
    data_seurat <- readRDS(argv$rds)
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
}

# integration --
anchor_features <- SelectIntegrationFeatures(object.list = data_seurat_list)
anchors <- FindIntegrationAnchors(object.list = data_seurat_list, anchor.features = anchor_features)
data_integrate <- IntegrateData(anchorset = anchors)
DefaultAssay(data_integrate) <- "integrated"

# clustering --
data_integrate <- data_integrate %>% 
  ScaleData() %>%
  RunPCA(verbose =F) %>%
  FindNeighbors(dims = 1:npcs, reduction = "pca", verbose = F) %>%
  FindClusters(resolution = resolution,verbose = F) %>%
  RunTSNE(dims=1:npcs, do.fast = TRUE, check_duplicates = FALSE, verbose = F) %>%
  RunUMAP(dims=1:npcs, reduction = "pca", verbose = F)

saveRDS(data_integrate, paste0(outdir, "/data_integrate.rds"))