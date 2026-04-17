
suppressWarnings(suppressMessages({
    library(argparser)
    library(Seurat)
    library(Matrix)
    library(tidyverse)
}))

dirnameScript <- function(){
    # get full directory path of current script located
    cmd = commandArgs(trailingOnly = FALSE)
    scriptName = sub('--file=', "", cmd[grep('^--file=', cmd)])
    if (length(scriptName) > 0) {
        path = normalizePath(scriptName)
        dirname = dirname(path)
    } else {
        print('Warning: not a runtime environment, using current directory instead.')
        dirname = getwd()
    }
    return(dirname)
}
 
sourceScript <- function(x, dir='/SGRNJ06/randd/USER/wangjingshen/script/') {
    # try source script from directory one by one: file itself, dir from argument, script directory
    scriptdir = dirnameScript()
    searchpath = c(x,
                   file.path(dir, x),
                   file.path(scriptdir, x))
    sourced = FALSE
    for (i in searchpath) {
        if (file.exists(i)) {source(i); print(i); sourced = TRUE; break}
    }
    if (!sourced) {stop(paste0('can not source: ', x))}
}
scriptdir = dirnameScript() #

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--matrix_10X", help = "mix matrix 10X")
#argv <- add_argument(argv, "--matrix_10X_gene", help = "mix matrix 10X gene")
argv <- add_argument(argv, "--percent_threshold", help = "percent threshold, Default: 0.9")
#argv <- add_argument(argv, "--outdir", help = "outdir, Default: outdir/")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

matrix_10X <- argv$matrix_10X
matrix_10X_gene <- str_glue("{matrix_10X}/features.tsv.gz")
percent_threshold <- ifelse(is.na(argv$percent_threshold), 0.9, as.numeric(argv$percent_threshold))
#outdir <- ifelse(is.na(argv$outdir), "outdir/", argv$outdir)
name <- argv$name

if(!dir.exists(name)){
    dir.create(name, recursive = TRUE)
    dir.create(str_glue("{name}/{name}_human"), recursive = TRUE)
    dir.create(str_glue("{name}/{name}_mouse"), recursive = TRUE)
}


data_count <- Read10X(matrix_10X)

# read gene list 
hs_genes <- read.table(str_glue("{scriptdir}/../data/hs_genes.tsv"))[,1]
mm_genes <- read.table(str_glue("{scriptdir}/../data/mmu_genes.tsv"))[,1]

#
mtx_hs_genes <- intersect(row.names(data_count), hs_genes)  # hs gene
mtx_mm_genes <- intersect(row.names(data_count), mm_genes)  # mm gene
mtx_others <- setdiff(row.names(data_count), c(mtx_hs_genes, mtx_mm_genes)) # others
print(paste0("hs genes: ", length(mtx_hs_genes)))
print(paste0("mm genes: ", length(mtx_mm_genes)))
print(paste0("other genes: ", length(mtx_others)))

# get UMI percent
hs_percent <- colSums(data_count[ mtx_hs_genes,])/ colSums(data_count)
mm_percent <- colSums(data_count[ mtx_mm_genes,])/ colSums(data_count)

hs_bc <- colnames(data_count)[ hs_percent > percent_threshold & mm_percent < 1 - percent_threshold]
mm_bc <- colnames(data_count)[ mm_percent > percent_threshold & hs_percent < 1 - percent_threshold]

bc_df <- data.frame(barcode = colnames(data_count),
                    species = colnames(data_count))

bc_df$species[ bc_df$barcode %in% hs_bc] ="human"
bc_df$species[ bc_df$barcode %in% mm_bc] ="mouse"
bc_df$species[ !bc_df$barcode %in% c(hs_bc, mm_bc)] ="doublets"
write.table(bc_df, str_glue("{name}/{name}_species.tsv"),sep="\t",quote=F,row.names=F)


function_subset_mtx <- function(species){
    if(species == "human"){
        subset_genes <- mtx_hs_genes
    }else{
        subset_genes <- mtx_mm_genes       
    }
    # mtx --
    data_subset <- data_count[subset_genes, colnames(data_count) %in% bc_df$barcode[ bc_df$species==species]]
    #data_subset <- data_subset[ rowSums(data_subset)>0,]    
    writeMM(obj = data_subset, file = str_glue("{name}/{name}_{species}/matrix.mtx"))
    
    # gene --
    genes_raw <- read.table(matrix_10X_gene, sep="\t")
    print(genes_raw[1:2,])
    genes_subset <- genes_raw[ match(row.names(data_subset), genes_raw$V2),]

    if(identical(genes_subset$V2, row.names(data_subset))){
        print("genes.tsv and matrix.mtx genes: True")
        write.table(genes_subset, file = str_glue("{name}/{name}_{species}/genes.tsv"), quote = F, sep = '\t', col.names = F, row.names = F)
    }else{
        print("genes.tsv and matrix.mtx genes: False, please check.")
        quit()
    }
    # barcode
    write.table(colnames(data_subset), file = str_glue("{name}/{name}_{species}/barcodes.tsv"), quote = F, col.names = F, row.names = F)
}
function_subset_mtx("human")
function_subset_mtx("mouse")


print(str_glue('{name} mix determination done.'))