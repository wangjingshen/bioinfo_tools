# Created on 2024.07.03
# @author: wangjingshen
# http://192.168.2.19:10003/notebooks/rd_project/invadeseq/data/RD20102301/16s/2024-07-03_make_res.ipynb
# http://192.168.2.19:10003/notebooks/rd_project/invadeseq/data/RD20102301/16s/2024_07_03_reproduce_genus_cell.ipynb

suppressWarnings(suppressMessages({
    library(argparser)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--samplename", help = "sample name")
argv <- add_argument(argv, "--rna_dir", help = "rna dir")
argv <- add_argument(argv, "--pathseq_dir", help = "pathseq dir")
argv <- parse_args(argv)

samplename <- argv$samplename
rna_dir <- argv$rna_dir
pathseq_dir <- argv$pathseq_dir
outdir <- paste0(argv$pathseq_dir, "/outdir/")

# mkdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

cat(paste0("UMI_df: sample id: ", samplename, "\n"))
cat("---\n")
## read file ----
# V1:tax_id; V2:taxonomy; ... 
df_16s_id_taxonomy <- read.table(paste0(pathseq_dir, "/", samplename, "_pathseq_score.txt"), sep="\t", header = T)
# V1: reads id; V2: bc; V3:umi; V4:id; V5:mapping quality; V6:genus;
df_16s_id_umi <- read.table(paste0(pathseq_dir, "/", samplename, ".invadeseq.raw.readnamepath"), sep="\t")
bc_rna <- read.table(paste0(rna_dir, "/outs/filtered/barcodes.tsv.gz"))[,1]  # get bc
df_16s_id_umi_rna <- df_16s_id_umi [ df_16s_id_umi$V2 %in% bc_rna,]

## process taxonomy
tax_id_detect <- unique(unlist(str_split(df_16s_id_umi_rna$V4,pattern = ",")))
df_16s_id_taxonomy <- df_16s_id_taxonomy[ df_16s_id_taxonomy$tax_id %in% tax_id_detect,]   # detect
taxonomy_detect <- str_split(df_16s_id_taxonomy$taxonomy,pattern = "\\|",simplify = F)
# df_16s_id_taxonomy$taxonomy[lapply(taxonomy_detect,function(x){sum(duplicated(x))}) > 0][1]
# root|cellular_organisms|Bacteria|Terrabacteria_group|Actinobacteria|Actinobacteria|Corynebacteriales|Nocardiaceae|Nocardia|Nocardia_sp._BMG51109
taxonomy_detect <- lapply(taxonomy_detect, unique)  # del repeated fields
# fill missing fields
taxonomy_detect[ unlist(lapply(taxonomy_detect, length)) == 8] <- lapply(taxonomy_detect[ unlist(lapply(taxonomy_detect, length))==8] , function(x){x[c(1:8,8)]}) 
#table(unlist(lapply(taxonomy_detect, length)))

taxonomy_detect <- sapply(taxonomy_detect, function(x){ paste0(x[c(1:9)], collapse = "|") })  # select top 9 fields
table(unlist(lapply(str_split(taxonomy_detect, pattern = "\\|", simplify = F),length)))
df_16s_id_taxonomy$taxonomy <- taxonomy_detect


function_add_taxonomy <- function(row_id){
    df_merge <- df_16s_id_taxonomy[ df_16s_id_taxonomy$tax_id %in% unlist(str_split(df_16s_id_umi_rna[row_id, 4], pattern = ",")), ]
    df_merge$barcode <- df_16s_id_umi_rna[row_id, 2]
    df_merge$umi <- df_16s_id_umi_rna[row_id, 3]
    df_merge$reads_id <- df_16s_id_umi_rna[row_id,1]
    df_merge$type <- "species"  # fix
    return(df_merge)
}
df_16s_id_taxonomy_umi_rna <- do.call(rbind, lapply(1:nrow(df_16s_id_umi_rna), function_add_taxonomy))


function_stat <- function(df){
    level <- unique(df$type)

    # output umi
    res_umi <- df %>% group_by(barcode, taxonomy, type) %>%
        summarise(umi_count = length(unique(umi)),
                  reads_count = length(unique(reads_id)))
    res_umi <- res_umi[ order(res_umi$umi_count, res_umi$reads_count, decreasing = c(T,T)), ]
    write.table(res_umi, paste0(outdir, "/", samplename, "_pathseq_", level, "_umi.tsv"), sep = "\t", quote = F, row.names = F)

    # output reads
    res_reads <- df %>% group_by(barcode, umi, taxonomy, type) %>%
        summarise(reads_count = length(unique(reads_id)))
    res_reads <- res_reads[ order(res_reads$reads_count, decreasing = T), ]
    #write.table(res_reads, paste0(outdir, "/", samplename, "_pathseq_", level, "_reads.tsv"), sep = "\t", quote = F, row.names = F)

    # output bc
    res_bc <- df %>% group_by(taxonomy, type) %>%
        summarise(barcode_count = length(unique(barcode)))
    res_bc <- res_bc[ order(res_bc$barcode_count, decreasing = T), ]
    write.table(res_bc, paste0(outdir, "/", samplename, "_pathseq_", level, "_barcode_count.tsv"), sep = "\t", quote = F, row.names = F)

    return(list(res_umi = res_umi,
                res_reads = res_reads,
                res_bc = res_bc))
}
res_species <- function_stat(df_16s_id_taxonomy_umi_rna)
# output stat
stat_df <- data.frame(term = c("total_barcode","total_species","total_umi"),
                      count = c(length(unique(res_species$res_umi$barcode)),
                                length(unique(res_species$res_umi$taxonomy)),
                                sum(res_species$res_umi$umi_count)))
write.table(stat_df, paste0(outdir, "/", samplename, "_pathseq_species_stat.tsv"), sep = "\t", quote = F, row.names = F)


## level switch --
df_level_switch <- data.frame(level = c("kingdom","phylum","class","order","family","genus","species"),
                              n = c(3:9))
function_switch_level <- function(df, level){
    df$taxonomy <- word(df$taxonomy, start = 1,end = df_level_switch[ df_level_switch$level == level, 2],sep = "\\|")
    df$type <- df_level_switch[ df_level_switch$level == level, 1]
    return(df)
}
df_genus <- function_switch_level(df_16s_id_taxonomy_umi_rna, "genus")
res_species <- function_stat(df_genus)


#cat(paste0("total barcode: ", length(unique(res_umi$barcode)), "\n"))
#cat(paste0("total taxonomy: ", length(unique(res_umi$taxonomy)), "\n"))
#cat(paste0("total umi: ", sum(res_umi$umi_count), "\n"))
#cat("Done.\n")
#cat("------\n")