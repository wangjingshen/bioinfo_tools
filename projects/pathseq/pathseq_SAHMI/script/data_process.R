suppressMessages(suppressWarnings({
    library(argparser)
    library(Seurat)
    library(tidyverse)
}))


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rna_path", help = "rna path")
argv <- add_argument(argv, "--fj_path", help = "fj path")
argv <- add_argument(argv, "--rna_spname", help = "rna_spname")
argv <- add_argument(argv, "--fj_spname", help = "fj_spname")
argv <- parse_args(argv)

outdir <- "./SAHMI_data/"

rna_path <- argv$rna_path
rna_spname <- argv$rna_spname
rna_mtx <- paste0(rna_path, "/", rna_spname, "/outs/filtered/") 
print(rna_mtx)

fj_path <- argv$fj_path
fj_spname <- argv$fj_spname
fj_mtx <- paste0(fj_path, "/SAHMI_result_RNA/6.taxa_count/", fj_spname, ".counts.txt")
bc_map <- paste0(fj_path, "/convert_data/", fj_spname, "/02.convert/barcode_convert.json")

# outdir <- ifelse(is.na(argv$outdir), "./16S_stat/", argv$outdir)
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

## read data
data_RNA <- Read10X(rna_mtx)
data_16S <- read.table(fj_mtx, sep="\t", header = T)
colnames(data_RNA) <- gsub("_", "", colnames(data_RNA))

## read barcode map 
barcode_map <- read.table(bc_map, sep = "\t")
barcode_map <- barcode_map[ 2: nrow(barcode_map) - 1, , drop=F]
barcode_map <- gsub(" ", "", barcode_map[,1])
barcode_map <- gsub(",", "", barcode_map)
barcode_map <- as.data.frame(str_split(barcode_map, pattern = ":", simplify = T))
print(paste0("n cells of RNA: ", ncol(data_RNA)))
print(paste0("n cells of intersection: ", length(intersect(barcode_map$V2, colnames(data_RNA)))))
if(ncol(data_RNA) != length(intersect(barcode_map$V2, colnames(data_RNA)))){
    print("Cells not in fj:")
    print(setdiff(colnames(data_RNA),barcode_map$V2))   # print self bc
    print(barcode_map[ match(setdiff(colnames(data_RNA),barcode_map$V2), barcode_map$V2),1])  # print 10X bc; not in bc 
}

barcode_map_bc <- barcode_map[ barcode_map$V2 %in% colnames(data_RNA),]
data_16S_bc <- data_16S[ data_16S$barcode %in% barcode_map_bc$V1,]
barcode_map_bc_16S <- barcode_map_bc[ barcode_map_bc$V1 %in% data_16S_bc$barcode,]

## trans 10X bc to sgr
data_16S_bc$barcode_10X <- data_16S_bc$barcode
data_16S_bc$barcode <- plyr::mapvalues(data_16S_bc$barcode_10X, from = barcode_map_bc_16S$V1, to = barcode_map_bc_16S$V2)
write.table(data_16S_bc, paste0(outdir, "/", fj_spname, "_counts_bc.tsv"), sep = "\t", quote = F, row.names = F)

data_16S_bc[is.na(data_16S_bc)] <- 0  # NA to 0

function_nunique <- function(x){
    x <- x[x!='0']  # del 0(NA)
    return(unique(x))
}

## stat genus
stat_df_genus <- data_16S_bc %>% 
    group_by(genus) %>%
    filter(genus!='0') %>%
    filter(rank=='g') %>%
    summarise(detect_bc = length(function_nunique(barcode)),
              total_bc = ncol(data_RNA),
              mean_umi_in_detect_bc = sum(counts)/sum(counts>0)) %>% 
    arrange(desc(detect_bc)) %>%
    mutate(detect_percent = paste0(round(100*detect_bc/total_bc,2), " %"))

# output
#write.table(stat_df_detected, paste0(outdir, "/", fj_spname, "_detected.tsv"), sep = "\t", quote = F, row.names = F)
#write.table(stat_df_species, paste0(outdir, "/", fj_spname, "_species.tsv"), sep = "\t", quote = F, row.names = F)
write.table(stat_df_genus, str_glue("{outdir}/{fj_spname}_genus_detect_stat.tsv"), sep = "\t", quote = F, row.names = F)


# genus tsv
df_genus <- data_16S_bc[ data_16S_bc$rank =="g",]

df_genus <- spread(df_genus[,c("barcode","genus","counts")], key = "genus", value = "counts")
df_genus[is.na(df_genus)] <- 0

df_genus <- data.frame(barcode = colnames(data_RNA)) %>% 
    left_join(df_genus, by = c("barcode" = "barcode"))
if(!identical(df_genus$barcode, colnames(data_RNA))){
    print("bc not equal.")
    print(df_genus$barcode[1:5])
    print(colnames(data_RNA)[1:5])
}else{
    print("bc equal.")
}

df_genus[is.na(df_genus)] <- 0
write.table(df_genus, str_glue("{outdir}/{fj_spname}_genus_umi.tsv"), sep="\t", row.names = F, quote = F)

bc_sub_len <- nchar(df_genus$barcode[1])/3
print(paste("bc sub length:", bc_sub_len))

function_split_bc <- function(x){
    return(paste(rna_spname,
                 substring(x, 1,bc_sub_len),
                 substring(x, 1+bc_sub_len, 2*bc_sub_len),
                 substring(x, 1+2*bc_sub_len, 3*bc_sub_len), sep="_"))
}
df_genus$barcode <- sapply(df_genus$barcode, function_split_bc)
write.table(df_genus, str_glue("{outdir}/{fj_spname}_genus_umi_for_down_analysis.tsv"), sep="\t", row.names = F, quote = F)
