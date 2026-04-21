
options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--zl_matrix", help = "zl matrix")
argv <- add_argument(argv, "--virus_matrix", help = "virus matrix")
argv <- add_argument(argv, "--outdir", help = "outdir to save fig. Default: . ")
argv <- parse_args(argv)

if(is.na(argv$outdir)){
    argv$outdir = "."
}
if(!is.na(argv$rds)){   # argv$rds != "None", for python run R
    print(paste0("rds: ", argv$rds))
    data_seurat <- readRDS(argv$rds)
}
if(!is.na(argv$zl_matrix)){
    cat(paste0("zl matrix : ", argv$zl_matrix))
    data_seurat <- CreateSeuratObject(counts = Read10X(argv$zl_matrix))
    data_seurat$percent.mt <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
    data_seurat <- data_seurat %>% 
        NormalizeData() %>%
        FindVariableFeatures(nfeatures = 2000) %>%
        ScaleData()%>% 
        RunPCA(verbose =F) %>%
        FindNeighbors(dims = 1:20,verbose = F) %>%
        FindClusters(resolution = 1,verbose = F) %>%
        RunTSNE(dims=1:20, do.fast = TRUE, check_duplicates = FALSE, verbose = F) %>%
        RunUMAP(dims=1:20,verbose = F) 
}
cat(paste0("virus matrix : ", argv$virus_matrix))
virus_matrix <- as.data.frame(t(as.matrix(Read10X(argv$virus_matrix))))
print( iden )
data_seurat$EBER1 <- virus_matrix$EBER1
data_seurat$EBER2 <- virus_matrix$EBER2
data_seurat$EBNA2 <- virus_matrix$EBNA2
data_seurat$ZEBRA <- virus_matrix$ZEBRA
data_seurat$EBNA1 <- virus_matrix$EBNA1

# Calculate QC variables ----
p1 <- FeaturePlot(data_seurat, c("EBER1","EBER2","EBNA2","EBNA1","ZEBRA"),ncol = 3)
ggsave(paste0(argv$outdir, "/EBV_gene_plot.png"), plot = p1, height = 7, width = 12)

cat("----------------\nEBV plot done.")