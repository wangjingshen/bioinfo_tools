suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--matrix_10X", help = "matrix_10X")
argv <- add_argument(argv, "--fj_UMI", help = "fj UMI")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- add_argument(argv, "--name", help = "sample name, code will add _ZL and _FJ")
argv <- parse_args(argv)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = T)
}

function_mtx2seurat <- function(matrix_10X){
    # clustering --
    data_seurat <- CreateSeuratObject(Read10X(matrix_10X)) %>% 
      NormalizeData() %>%
      FindVariableFeatures(nfeatures = 2000) %>%
      ScaleData()%>% 
      RunPCA(verbose =F) %>%
      FindNeighbors(dims = 1:20,verbose = F) %>%
      FindClusters(resolution = 1,verbose = F) %>%
      RunUMAP(dims=1:20,verbose = F)
    return(data_seurat)
}

data_seurat <- function_mtx2seurat(argv$matrix_10X)
fj_UMI <- read.table(argv$fj_UMI, sep=",", header = T)
zl_gene_list <- intersect(row.names(data_seurat), c("ADAMTS6","ANKRD33B","CCND1","CD247","GRK3","HMOX1","MTA1","RXRA","TBC1D8","SOX11"))
fj_gene_list <- intersect(colnames(fj_UMI), c("ADAMTS6","ANKRD33B","CCND1","CD247","GRK3","HMOX1","MTA1","RXRA","TBC1D8","SOX11"))

for( i in fj_gene_list){
    data_seurat@meta.data[,paste0("FJ_", i)] <- fj_UMI[,i]
}

# plot ----
function_gene_plot <- function(data, gene, name){
    FeaturePlot(data, gene) + ggtitle(paste0(name, ": ", gene))
    ggsave(paste0(argv$outdir, "/", name, "_", gene, ".png"),height = 5,width = 6) 
}
sapply(zl_gene_list, function_gene_plot, data = data_seurat, name = paste0(argv$name, "_ZL"))

fj_gene_list_update <- paste0("FJ_", fj_gene_list)
sapply(fj_gene_list_update, function_gene_plot, data = data_seurat, name = argv$name )


# stat df ----
function_stat_df <- function(data, gene_list, mod, name){
    if(mod == "zl"){
        res <- data.frame(genes = gene_list,
                          express_cells = sapply(gene_list, function(x){sum(data@assays$RNA@counts[x,]>0)}),
                          total_cells = rep(ncol(data), length(gene_list))) %>%
        mutate(express_percent = paste0(round(express_cells/total_cells,4)*100, " %"))
        write.table(res, paste0(argv$outdir, "/", name ,"_ZL_stat_df.tsv"), sep="\t", quote = F, row.names = F) 
    }else{
        res <- data.frame(genes = gene_list,
                          express_cells = sapply(gene_list, function(x){sum(data@meta.data[,x]>0)}),
                          total_cells = rep(ncol(data), length(gene_list))) %>%
        mutate(express_percent = paste0(round(express_cells/total_cells,4)*100, " %"))
        write.table(res, paste0(argv$outdir, "/", name ,"_FJ_stat_df.tsv"), sep="\t", quote = F, row.names = F) 
    }
}

function_stat_df(data_seurat, zl_gene_list, 'zl', argv$name)
function_stat_df(data_seurat, fj_gene_list_update, 'fj', argv$name)
