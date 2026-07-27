suppressWarnings(suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
}))

color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--anno", help = "anno")
argv <- add_argument(argv, "--subset", help = "subset or not")
argv <- add_argument(argv, "--df_genus", help = "df genus")
argv <- add_argument(argv, "--rna_spname", help = "rna sample name")
argv <- add_argument(argv, "--spname", help = "sample name")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

rds <- argv$rds
df_genus <- argv$df_genus
rna_spname <- argv$rna_spname
spname <- argv$spname
outdir <- argv$outdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}


# read RDS -- 
data_seurat <- readRDS(rds)

if(!is.na(argv$subset)){
    data_seurat <- subset(data_seurat, subset= sample == argv$subset)
}
if(!is.na(argv$anno)){
    data_seurat$cluster <- data_seurat$cluster_singleR
}
data_seurat$spname <- spname

# 16S
df_genus <- read.table(df_genus, sep="\t", header = T)

# check
if( !identical(colnames(data_seurat), df_genus$barcode) ){
    print("barcode not equal.")
    print(colnames(data_seurat)[1:5])
    print(df_genus$barcode[1:5])
    quit()
}else{
    print("barcode equal.")
}

data_seurat$total_bacteria <- rowSums(df_genus[,-1])

# plot
options(repr.plot.height = 6, repr.plot.width = 7)
FeaturePlot(data_seurat, features = "total_bacteria", pt.size = 1,order = T)
ggsave(str_glue('{outdir}/Featureplot_total_bacteria_{spname}.png'), width = 7, height = 5)


stat_bc <- data_seurat@meta.data %>% 
    group_by(spname, cluster) %>%
    summarise(detect_bc = sum(total_bacteria > 0),
              total_bc = length(total_bacteria),
              mean_umi_in_detect_bc = sum(total_bacteria)/sum(total_bacteria>0) ) %>%
    mutate(detect_percent = paste0(round(100*detect_bc/total_bc,2), " %"))

write.table(stat_bc, str_glue("{outdir}/bacteria_detect_stat_{spname}.tsv"),sep="\t", row.names = F, quote = F)

print('Done')