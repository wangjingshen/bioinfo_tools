#! /usr/bin/env Rscript
suppressWarnings(suppressMessages({
    library(argparser)
    library(sceasy)
    library(Seurat)
    library(scater)
    library(patchwork)
    library(reticulate)
    library(tidyverse)
    use_condaenv('r4.1_env')
    loompy <- reticulate::import('loompy')
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--input", help = "input")
argv <- add_argument(argv, "--from_type", help = "from_type, anndata or seurat or sce or loom")
argv <- add_argument(argv, "--to_type", help = "to_type, anndata or loom or sce or seurat or cds")
argv <- add_argument(argv, "--outdir", help = "outdir, Default: ./outdir/")
argv <- parse_args(argv)


input <- argv$input
from_type <- argv$from_type
to_type <- argv$to_type
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

## trans to anndata --
# sce,loom,seurat -> anndata
if(to_type == "anndata"){
    if(from_type %in% c("seurat", "sce")){
        input <- readRDS(input)
    }
    convertFormat(input, from = from_type, to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.h5ad'))
}

## trans to seurat --
# 1) sce,loom -> anndata -> seurat; 2) others -> seurat
if(to_type == "seurat"){
    if(from_type %in% c("sce","loom")){
        if(from_type == "sce"){
            input <- readRDS(input)
        }
        convertFormat(input, from = from_type, to = "anndata", outFile = str_glue('{outdir}/tmp.h5ad'))
        data <- convertFormat(str_glue('{outdir}/tmp.h5ad'), from = "anndata", to = to_type)
        file.remove(str_glue('{outdir}/tmp.h5ad'))  # rm tmp 
    }else{
        data <- convertFormat(input, from = "anndata", to = to_type, main_layer = 'counts')
    }
    data@assays$RNA@key <- "rna_"  # fix
    saveRDS(data, file = str_glue('{outdir}/{from_type}_to_{to_type}.rds'))
}

## trans to loom --
# 1) seurat -> sce -> loom; 2) anndata -> seurat -> sce -> loom; 3) others -> loom
if(to_type == "loom"){
    if( from_type == "seurat"){
        input <- convertFormat(readRDS(input), from = from_type, to = "sce")
        convertFormat(input, from = "sce", to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.loom'))
    }
    if( from_type == "anndata"){
        input <- convertFormat(convertFormat(input, from = from_type, to = "seurat"), from = "seurat", to = "sce")
        convertFormat(input, from = "sce", to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.loom'))
    }
    if( from_type == "sce"){
        convertFormat(readRDS(input), from = from_type, to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.loom'))
    }
}

## trans to sce --
# 1) anndata -> seurat -> sce; 2) others -> sce
if(to_type == "sce"){
    if(from_type == "anndata"){
        input <- convertFormat(input, from = from_type, to = "seurat")
        convertFormat(input, from = "seurat", to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.rds'))
    }else{
        if(from_type == "seurat"){
            input <- readRDS(input)
        }
        convertFormat(input, from = from_type, to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.rds'))
    }
}

## trans to cds --
# sce,loom,seurat -> anndata -> cds
if(to_type == "cds"){
    if(from_type %in% c("seurat", "sce")){
        input <- readRDS(input)
    }
    convertFormat(input, from = from_type, to = "anndata", outFile = str_glue('{outdir}/tmp.h5ad'))
    convertFormat(str_glue('{outdir}/tmp.h5ad'), from = "anndata", to = to_type, outFile = str_glue('{outdir}/{from_type}_to_{to_type}.rds'))
    if(from_type %in% c("seurat", "sce")){
        file.remove(str_glue('{outdir}/tmp.h5ad'))  # rm tmp 
    }
}


print('convert finished')