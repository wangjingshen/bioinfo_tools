# down analysis for /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16s_sgr/script/pipeline.py
# total genus
# top10 genus
# cluster genus mean umi

suppressMessages(suppressWarnings({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(patchwork)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--path", help = "path")
argv <- add_argument(argv, "--spname", help = "spname")
argv <- add_argument(argv, "--outdir", help = "output dir, default: outdir")
argv <- parse_args(argv)

path <- argv$path
spname <- unlist(str_split(argv$spname, ','))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

function_cal <- function(spname){
    df_genus <- data.table::fread(str_glue("{path}/{spname}/outs/{spname}_raw_UMI_matrix.tsv.gz"))

    df <- data.frame(
        genus = df_genus[,1],
        positive_cells = rowSums(df_genus[, -1]>0),
        total_cells = ncol(df_genus)-1) %>%
        mutate(positive_percent = paste0(round(100*positive_cells/total_cells, 2)," %"))
    colnames(df)[1] <- "genus"

    df <- df[ order(df$positive_cells, decreasing=T),]
    write.table(df, str_glue("{outdir}/{spname}_stat.tsv"), sep="\t", quote = F, row.names = F)
}
lapply(spname, function_cal)

print('Done')