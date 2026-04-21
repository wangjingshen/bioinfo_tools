# Created on 2024.07.03
# @author: wangjingshen
# http://192.168.2.19:10003/notebooks/rd_project/invadeseq/data/RD20102301/16s/2024-07-03_make_res.ipynb
# http://192.168.2.19:10003/notebooks/rd_project/invadeseq/data/RD20102301/16s/2024_07_03_reproduce_genus_cell.ipynb
# taxid to tax

# init values
PIPE_topn = 10
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# import packages
suppressWarnings(suppressMessages({
    library(argparser)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--sample", help = "sample name")
argv <- add_argument(argv, "--match_dir", help = "match rna dir")
argv <- add_argument(argv, "--pathseq_dir", help = "pathseq dir")
argv <- add_argument(argv, "--outdir", help = "outdir")   #04
argv <- parse_args(argv)

sample <- argv$sample
match_dir <- argv$match_dir
pathseq_dir <- argv$pathseq_dir
outdir <- argv$outdir

cat(paste0("genus_df: sample id: ", sample, "\n"))
# mkdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

##
data_filter <- read.table(str_glue('{pathseq_dir}/{sample}.invadeseq.genus.csv'), sep=",", header = T)
match_bc <- read.table(str_glue('{match_dir}/outs/filtered/barcodes.tsv.gz'))
match_bc <- paste0(sample, "_", match_bc[,1])
data_filter <- data.frame(barcode = match_bc) %>%
        left_join(data_filter , by = c("barcode" = "barcode"))
data_filter[is.na(data_filter)] <- 0
write.table(data_filter, str_glue('{outdir}/{sample}_genus_umi.tsv'), sep = "\t", quote = F, row.names = F)

# genus barcode count
genus_bc <- data.frame(genus = names(sort(colSums(data_filter[,-1] > 0), decreasing = T)),    # -1 rm barcode
                       detect_barcode = as.numeric(sort(colSums(data_filter[,-1] > 0), decreasing = T))) %>%
                       mutate(cell_number = length(match_bc),
                              percent = paste(round(detect_barcode/cell_number*100, 2),"%"))
write.table(genus_bc, str_glue('{outdir}/{sample}_genus_detect_barcode.tsv'), sep = "\t", quote = F, row.names = F)


# genus umi
genus_umi <- data.frame(genus = names(colSums(data_filter[,-1])),    # -1 rm barcode
                        umi = as.numeric(colSums(data_filter[,-1]))) %>%
    arrange(desc(umi))  %>%
    mutate(total_umi = sum(umi),
           percent = paste(round(umi/total_umi*100, 2),"%"))
write.table(genus_umi, str_glue('{outdir}/{sample}_genus_detect_umi.tsv'), sep = "\t", quote = F, row.names = F)

cat("Done.\n")
cat("------\n")