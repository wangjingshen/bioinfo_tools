# human
suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args ----
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--version", help = "version")
argv <- add_argument(argv, "--spname", help = "spname")
argv <- add_argument(argv, "--tag_file", help = "tag file")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: outdir")
argv <- parse_args(argv)

rds <- argv$rds
version <- argv$version
spname <- unlist(strsplit(argv$spname, split = ","))
tag_file <- unlist(strsplit(argv$tag_file, split = ",")) 
outdir <- argv$outdir
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)
}

data_seurat <- readRDS(rds)
if(version == "singleR"){
    data_seurat$cluster <- data_seurat$celltype
}

# read tag
function_read_tag_csv <- function(i){
    tag <- read.table(tag_file[i], sep="\t", header = T, row.names = 1, check.names = F)
    if(version != "singleR"){
        row.names(tag) <- paste0(spname[i], "_", row.names(tag))
    }
    return(tag[,-c(1:4)]) # del tSNE_1  tSNE_2  cluster Gene_Counts 
}
tag <- do.call(rbind, lapply(1:length(spname), function_read_tag_csv))
tag <- tag[ colnames(data_seurat),]

print(dim(data_seurat))
print(dim(tag))

data_seurat@meta.data <- cbind.data.frame(data_seurat@meta.data, tag)


options(repr.plot.height = 6, repr.plot.width = 9)
DimPlot(data_seurat, group.by = "cluster", label = T, repel = T, cols = color_protocol)
ggsave(str_glue("{outdir}/Featureplot_cluster.png"), width = 7, height = 6)


for (plot_tag in colnames(tag)[ -ncol(tag)]){
    FeaturePlot(data_seurat, features = plot_tag, pt.size=1, order = T)
    ggsave(str_glue("{outdir}/Featureplot_{plot_tag}.png"), width = 7, height = 5)
    
    data_seurat$analysis_col <- data_seurat@meta.data[, plot_tag]
    df_stat <- data_seurat@meta.data %>% 
        group_by(cluster) %>%
        summarise(tag = plot_tag,
                  ncell_total = length(analysis_col),
                  ncell_exp = sum(analysis_col > 0),
                  umi_total = sum(analysis_col),
                  umi_mean = round(mean(analysis_col),2)) %>%
        mutate(ncell_exp_pct = paste0(round(100 * ncell_exp/ncell_total, 2), " %"))
    #df_stat[is.na(df_stat)] <- 0
    write.table(df_stat[,c(1:4,7,5:6)], str_glue("{outdir}/detect_stat_{plot_tag}.tsv"), sep="\t", row.names = F, quote = F)
}


df_mean <- do.call(cbind,lapply(colnames(tag)[ -ncol(tag)],function(x){
    data_seurat$analysis_col <- data_seurat@meta.data[, x]
    df_mean <- data_seurat@meta.data %>% 
        group_by(cluster) %>%
        summarise(mean = round(mean(analysis_col),2))
    x < gsub("-","_",x)
    colnames(df_mean)[2] <- str_glue("mean_{x}")
    return(df_mean)
}))

write.table(df_mean, str_glue("{outdir}/detect_stat_mean.tsv"), sep="\t", row.names = F, quote = F)