# Created on 2024.07.03
# @author: wangjingshen

# init values
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
primer_file <- "/SGRNJ06/randd/PROJECT/RD20102301_DZH/16S_display/primer/fasta_file/16S_primer_degenerate.fasta"

# import packages
suppressWarnings(suppressMessages({
    library(argparser)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--df", help = "df")
argv <- add_argument(argv, "--name", help = "name")
argv <- add_argument(argv, "--top_n", help = "top_n")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- parse_args(argv)

df <- argv$df
name <- argv$name
top_n <- argv$top_n
outdir <- argv$outdir

# mkdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}


## read file ----
df <- read.table(df, sep="\t", header = 0)
colnames(df) <- "primer"
df$count = 1
df_stat <- df %>% 
    group_by(primer) %>%
    summarise(primer_count = sum(count)) %>%
    mutate(total_count = sum(primer_count),
           percent = paste0(round( 100 * primer_count/total_count, 2), " %")) %>%
    arrange(desc(primer_count))
    
write.table(df_stat, str_glue('{outdir}/{name}_top{top_n}bp_stat.tsv'), sep = "\t", quote = F, row.names = F)

# update
function_update <- function(p_update){
    p_update <- paste0("^",p_update)
    df_update <- df_stat[grep(p_update, df_stat$primer),]
    stat_update <- data.frame(primer = gsub("\\^","",p_update),
                              primer_count = sum(as.numeric(df_update$primer_count)),
                              total_count = as.numeric(df_update$total_count[1]),
                              percent = paste0(round(100*sum(df_update$primer_count)/ df_update$total_count[1],2), " %"))
    return(stat_update)
}


primer <- read.table(primer_file, sep="\t")
primer <- primer[ -grep(">", primer$V1),]
df_primer <- do.call(rbind, lapply(primer, function_update))
df_primer <- df_primer[ order(df_primer$primer_count, decreasing = T),]
    
write.table(df_primer, str_glue('{outdir}/{name}_primer_stat.tsv'), sep = "\t", quote = F, row.names = F)

cat("Done.\n")
cat("------\n")