suppressWarnings(suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(patchwork)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--cele_path", help = "cele_path")
argv <- add_argument(argv, "--species", help = "species, human or mouse. Default: human")
argv <- add_argument(argv, "--name", help = "name")
argv <- add_argument(argv, "--topn", help = "topn. Default: 20")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: stat/")
argv <- parse_args(argv)

cele_path = argv$cele_path
species <- ifelse(is.na(argv$species), "human", argv$species)
if(species == "human"){
    gene_trans <- "/SGRNJ06/randd/USER/wangjingshen/share/genome/celescope_v2/hs/homo.gtf.trans"
}
if(species == "mouse"){
    gene_trans <- "/SGRNJ06/randd/USER/wangjingshen/share/genome/celescope_v2/mm/Mus.gtf.trans"
}

name = argv$name

topn <- ifelse(is.na(argv$topn), 20, as.numeric(argv$topn))
outdir <- ifelse(is.na(argv$outdir), "stat/", argv$outdir)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)
}

# reads stat --
function_reads_stat <- function(mtx, name){
    reads_mtx <- read.table(mtx, sep="\t", header = T, row.names = 1)
    colnames(reads_mtx)[6] <- "reads"
    reads_mtx$percent <- reads_mtx$reads / sum(reads_mtx$reads)
    reads_stat <- reads_mtx[order(reads_mtx$percent, decreasing = T)[1:topn],c("reads","percent")]
    
    gene_trans <- read.table(gene_trans, sep=" ")[,c(2,4)]  # 2 ENS 4 gene
    reads_stat$gene <- gene_trans[match(row.names(reads_stat), gene_trans$V2), 2]
    reads_stat$percent <- paste0(round(100*reads_stat$percent, 4), " %")
    write.table(reads_stat[,c(3,1,2)], file = paste0(outdir, "/reads_stat_", name, ".xls"), sep="\t", quote = F, row.names = F) # gene,reads,percent
    #return(reads_stat[,c(3,1,2)])
    #return(reads_mtx)
}
function_reads_stat(mtx = name, name = name)


# umi stat
function_umi_stat <- function(mtx, name){
    umi_mtx <- Read10X(mtx)
    umi_stat <- data.frame(gene = row.names(umi_mtx),
                           UMI = rowSums(umi_mtx),
                           percent = rowSums(umi_mtx)/sum(rowSums(umi_mtx)))
    umi_stat <- umi_stat[order(umi_stat$percent, decreasing = T),]
    umi_stat$percent <- paste0(round(100*umi_stat$percent, 4), " %")
    write.table(umi_stat[1:topn,], file = paste0(outdir, "/UMI_stat_", name, ".xls"), sep="\t", quote = F, row.names = F)
    #return(umi_stat)
}
function_umi_stat(mtx = paste0(cele_path, "/", name, "/outs/raw/"), name = paste0(name, "_raw"))
function_umi_stat(mtx = paste0(cele_path, "/", name, "/outs/filtered/"), name = paste0(name, "_filtered"))

# end --
cat("stat done.\n")
cat("--------\n")