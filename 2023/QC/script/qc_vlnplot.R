
options(sci.pen = 10)
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
argv <- add_argument(argv, "--matrix_10X", help = "matrix, just for one sample")
argv <- add_argument(argv, "--species", help = "species, huamn or mouse or homo_mus")
argv <- add_argument(argv, "--sname", help = "sample name")
argv <- add_argument(argv, "--outdir", help = "outdir to save fig. Default: . ")
argv <- parse_args(argv)

if(!is.na(argv$rds)){   # argv$rds != "None", for python run R
    print(paste0("rds: ", argv$rds))
    data_seurat <- readRDS(argv$rds)
    Idents(data_seurat) <- argv$sname
}
if(!is.na(argv$matrix_10X)){
    print(paste0("matrix 10X : ", argv$matrix_10X))
    data_seurat <- CreateSeuratObject(counts = Read10X(argv$matrix_10X))
    Idents(data_seurat) <- argv$sname
}
if(!is.na(argv$outdir)){
    argv$outdir = "."
}
argv$species <- ifelse(is.na(argv$species), "homo_mus", argv$species)

# Calculate QC variables ----
# reference
# rb genes --
# https://cloud.tencent.com/developer/article/1605889?from=article.detail.1820062
# https://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query=ribosome&submit=Quick%0D%0ASearch
# HB --
# https://www.jianshu.com/p/207e71fc85fa
# https://www.genenames.org/tools/search/#!/?query=hemoglobin&rows=20&start=0&filter=locus_group:%22Protein-coding%20gene%22

if(argv$species == "human"){
    mt_pattern <- "^MT-"
    hb_genes <- read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/human_hemoglobin_genes.tsv", sep="\t")[,1]
    rb_genes <- read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/human_ribosomal_proteins_genes.tsv", sep="\t")[,1]
}
if(argv$species == "mouse"){
    mt_pattern <- "^Mt-"
    hb_genes <- read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/mouse_hemoglobin_genes.tsv", sep="\t")[,1]
    rb_genes <- read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/mouse_ribosomal_proteins_genes.tsv", sep="\t")[,1]
}
if(argv$species == "homo_mus"){
    mt_pattern <- "^M[Tt]-"
    hb_genes <- c(read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/human_hemoglobin_genes.tsv", sep="\t")[,1],
                  read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/mouse_hemoglobin_genes.tsv", sep="\t")[,1])
    rb_genes <- c(read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/human_ribosomal_proteins_genes.tsv", sep="\t")[,1],
                  read.table("/SGRNJ06/randd/USER/wangjingshen/script/QC/data/mouse_ribosomal_proteins_genes.tsv", sep="\t")[,1])
}

data_seurat$percent.mt <- PercentageFeatureSet(data_seurat, pattern = mt_pattern)
# rb genes --
rb.genes <- intersect(rownames(data_seurat@assays$RNA@counts), rb_genes)
#percent.rb <- Matrix::colSums(data_seurat@assays$RNA@counts[rb.genes,])/Matrix::colSums(data_seurat@assays$RNA@counts)*100
#data_seurat <- AddMetaData(data_seurat, percent.rb, col.name = "percent.rb")
data_seurat[["percent.rb"]]<-PercentageFeatureSet(data_seurat, features=rb.genes) 

# hb genes --
hb.genes <- intersect(rownames(data_seurat@assays$RNA@counts), hb_genes)
data_seurat[["percent.hb"]]<-PercentageFeatureSet(data_seurat, features=hb.genes) 

# plot ----
p1 <- VlnPlot(data_seurat, c("nFeature_RNA"),pt.size = 0) + theme(axis.title.x = element_blank(),legend.position = "none")
p2 <- VlnPlot(data_seurat, c("nCount_RNA"),pt.size = 0) + theme(axis.title.x = element_blank(),legend.position = "none")
p3 <- VlnPlot(data_seurat, c("percent.mt"),pt.size = 0) + theme(axis.title.x = element_blank(),legend.position = "none")
p4 <- VlnPlot(data_seurat, c("percent.rb"),pt.size = 0) + theme(axis.title.x = element_blank(),legend.position = "none")
p5 <- VlnPlot(data_seurat, c("percent.hb"),pt.size = 0) + theme(axis.title.x = element_blank(),legend.position = "none")

p <- patchwork::wrap_plots(list(p1,p2,p3,p4,p5),ncol = 5)
ggsave(paste0(argv$outdir, "/", argv$sname ,"_QC.png"), plot = p, height = 4.5, width = 10)

# print summary ----
print("nFeature_RNA summary:")
print(summary(data_seurat$nFeature_RNA))
print("nCount_RNA summary:")
print(summary(data_seurat$nCount_RNA))
print("percent.mt summary:")
print(summary(data_seurat$percent.mt))
print("percent.rb summary:")
print(summary(data_seurat$percent.rb))
print("percent.hb summary:")
print(summary(data_seurat$percent.hb))


cat("----------------\nQC done.")