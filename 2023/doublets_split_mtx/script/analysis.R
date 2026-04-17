
suppressWarnings(suppressMessages({
    library(argparser)
    library(Seurat)
    library(Matrix)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--mtx", help = "matrix 10X")
argv <- add_argument(argv, "--doublets_species", help = "doublets species")
#argv <- add_argument(argv, "--species1", help = "species1. Default: human")
#argv <- add_argument(argv, "--species2", help = "species2. Default: mouse")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- add_argument(argv, "--prefix", help = "prefix")
argv <- parse_args(argv)

#argv$species1 <- ifelse(is.na(argv$species1), "human", argv$species1)
#argv$species2 <- ifelse(is.na(argv$species2), "mouse", argv$species2)

if(is.na(argv$outdir)){
    argv$outdir = "./"
}
prefix <- argv$prefix
if(!dir.exists(paste0(argv$outdir, "/" ,argv$prefix, "_human"))){
    dir.create(paste0(argv$outdir, "/", argv$prefix, "_human"), recursive = TRUE)
}
if(!dir.exists(paste0(argv$outdir, "/" ,argv$prefix, "_mouse"))){
    dir.create(paste0(argv$outdir, "/" ,argv$prefix, "_mouse"), recursive = TRUE)
}

function_subset_mtx <- function(data_count, matrix_10X, doublets_species, outdir, prefix, species){
    # mtx --
    data_count_subset <- data_count[ ,colnames(data_count) %in% doublets_species$cell_id[ doublets_species$Species==species]]
    data_count_subset <- data_count_subset[ rowSums(data_count_subset)>0,]    
    #sum(rowSums(data_count_subset)==0)
    writeMM(obj = data_count_subset, file = paste0(outdir, "/", prefix, "_", species, "/matrix.mtx"))
    
    # gene name
    genes_raw <- read.table(paste0(matrix_10X, "/genes.tsv"),sep="\t")
    genes_subset <- genes_raw[ match(row.names(data_count_subset), genes_raw$V2),]
    #genes_subset <- genes_raw[ match(genes_raw$V2, row.names(data_count_subset)),]  # test, error
    #genes_subset <- genes_raw[ genes_raw$V2 %in% row.names(data_count_subset),]
    if(identical(genes_subset$V2, row.names(data_count_subset))){
        print("genes.tsv and matrix.mtx genes: True")
        write.table(genes_subset, file = paste0(outdir, "/", prefix, "_", species, '/genes.tsv'), quote = F, sep = '\t', col.names = F, row.names = F)
    }else{
        print("genes.tsv and matrix.mtx genes: False, please check.")
        quit()
    }

    # barcode
    write.table(colnames(data_count_subset), file = paste0(outdir, "/", prefix, "_", species, '/barcodes.tsv'), quote = F, col.names = F, row.names = F)
}

function_doublets_split_mtx <- function(matrix_10X, doublets_species, prefix, outdir){
    data_count <- Read10X(matrix_10X)
    doublets_species <- read.table(doublets_species,sep="\t",header = T)
    doublets_species[1,1] <- colnames(data_count)[1]

    if(identical(doublets_species$cell_id, colnames(data_count))){
        print("doublets cellid and mtx barcode: True")
        function_subset_mtx(data_count, matrix_10X, doublets_species, outdir, prefix, species = "human")
        function_subset_mtx(data_count, matrix_10X, doublets_species, outdir, prefix, species = "mouse")
    }else{
        print("doublets cellid and mtx barcode: False, please check.")
        quit()
    }
}

function_doublets_split_mtx(matrix_10X = argv$mtx,
                            doublets_species = argv$doublets_species,
                            prefix = prefix,
                            outdir = argv$outdir)

print('Done.')