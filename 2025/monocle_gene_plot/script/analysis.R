
suppressWarnings(suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(monocle)}
))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--cds", help = "cds, monocle cds")
argv <- add_argument(argv, "--mode", help = "mode, gene or pathway")
argv <- add_argument(argv, "--genes", help = "genes, split by ,")
argv <- add_argument(argv, "--df_ucell", help = "ucell tsv")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: HBV_stat")
argv <- parse_args(argv)

# default value --
#argv$gene <- ifelse(is.na(argv$gene), c("AB211369.1","BAE53654.1","BAE53655.1","BAE53656.1","BAE53658.1","BAE53661.1","BAE53662.1"), argv$gene)
mode <- unlist(strsplit(argv$mode, split = ","))
if(mode == "gene"){
    genes <- unlist(strsplit(argv$genes, split = ","))
    print(genes)
}
if(mode == "pathway"){
    df_ucell<- argv$df_ucell
    print(df_ucell)
}

outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

# read cds -- 
cds <- readRDS(argv$cds)

function_gene_plot <- function(cds, gene){
    cds$analysis_gene <- cds@assayData$exprs[gene,]
    if(sum(cds$analysis_gene > 0) >0){
        plot_cell_trajectory(cds, color_by = "analysis_gene") +
            scale_color_gradient(low = "lightgrey", high = "red") +
            labs(color = gene)
        ggsave(str_glue("{outdir}/monocle_gene_{gene}.png"), width = 6, height = 6)
    }else{
        plot_cell_trajectory(cds, color_by = "analysis_gene") +
            scale_color_gradient(low = "lightgrey", high = "lightgrey") +
            labs(color = gene)
        ggsave(str_glue("{outdir}/monocle_gene_{gene}.png"), width = 6, height = 6) 
    }
}

function_pathway_plot <- function(cds, df, pathway){
    cds$analysis_pathway <- df[,pathway]
    plot_cell_trajectory(cds, color_by = "analysis_pathway") +
        scale_color_gradient(low = "lightgrey", high = "red") +
        labs(color = pathway)
        pathway <- gsub(" ","_", pathway)
    ggsave(str_glue("{outdir}/monocle_pathway_{pathway}.png"), width = 6, height = 6)
}

if(mode == "gene"){
    out = lapply(genes, function_gene_plot, cds = cds)
}

if(mode == "pathway"){
    df_ucell <- as.data.frame(data.table::fread(df_ucell))
    row.names(df_ucell) <- df_ucell$V1
    df_ucell <- df_ucell[ row.names(cds@phenoData@data), ]

    if(identical(row.names(cds@phenoData@data), row.names(df_ucell))){
        pathways <- colnames(df_ucell)[(which(colnames(df_ucell)=="AB_status")+1) : ncol(df_ucell)]   # for
        out = lapply(pathways, function_pathway_plot, cds = cds, df = df_ucell)
    }else{
        print("bc not equal, please check.")
        print(row.names(cds@phenoData@data)[1:3])
        print(row.names(df_ucell)[1:3])
        quit()
    }
}

# end --
cat("Monocle plot done.\n")
cat("--------\n")