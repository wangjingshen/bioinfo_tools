# down analysis for /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16s_sgr/script/pipeline.py

suppressMessages(suppressWarnings({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(patchwork)
}))
top_n <- 5
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--anno", help = "anno")
argv <- add_argument(argv, "--subset", help = "subset")
argv <- add_argument(argv, "--df_genus", help = "df genus")
argv <- add_argument(argv, "--rna_spname", help = "rna spname name,split by ,")
argv <- add_argument(argv, "--spname", help = "spname name,split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

rds <- argv$rds
df_genus <- argv$df_genus
rna_spname <- argv$rna_spname
spname <- argv$spname
outdir <- argv$outdir

#
data_seurat <- readRDS(rds)
if(!is.na(argv$subset)){
    data_seurat <- subset(data_seurat, subset= sample == argv$subset)
}
if(!is.na(argv$anno)){
    data_seurat$cluster = data_seurat$cluster_singleR
}
# add barcode
data_seurat$spname <- spname

df_genus <- read.table(df_genus, sep="\t",header = T, row.names = 1)

# check
if( !identical(colnames(data_seurat), row.names(df_genus)) ){
    print("barcode not equal.")
    print(colnames(data_seurat)[1:5])
    print(row.names(df_genus)[1:5])
    quit()
}else{
    print("barcode equal.")
}

# plor cluster
options(repr.plot.width=7, repr.plot.height = 5)
DimPlot(data_seurat, group.by = "cluster")
ggsave(str_glue(outdir, "/cluster_",spname, ".png"), width = 7, height = 5)


data_seurat@meta.data <- cbind.data.frame(data_seurat@meta.data, df_genus)

top10_genus <- names(sort(colSums(df_genus > 0), decreasing = T)[1:10])


# all
options(repr.plot.width=16, repr.plot.height = 10)
FeaturePlot(data_seurat, names(sort(colSums(df_genus > 0), decreasing = T)[1:10]), pt.size=1, order = T)
ggsave(str_glue("{outdir}/{spname}_top10_genus.png"), width=16, height = 10)

# split
seq= 1
for (i in top10_genus){
    FeaturePlot(data_seurat, features = i, pt.size=1, order = T)
    ggsave(str_glue("{outdir}/{spname}_{seq}_{i}.png"), width=7, height = 5)
    
    data_seurat$analysis_col <- data_seurat@meta.data[,i]
    df_stat <- data_seurat@meta.data %>% 
        group_by(cluster) %>%
        summarise(genus = i,
                  detect_bc = sum(analysis_col >0),
                  total_bc = length(analysis_col),
                  mean_umi_in_detect_bc = sum(analysis_col)/sum(analysis_col>0) ) %>%
        mutate(detect_percent = paste0(round(100*detect_bc/total_bc,2), " %"))  
    df_stat[is.na(df_stat)] <- 0
    write.table(df_stat, str_glue("{outdir}/{spname}_{seq}_{i}_stat.tsv"),sep="\t",row.names=F,quote = F)
    seq = seq + 1
}
