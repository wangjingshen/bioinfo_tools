### v2 for lims
# fj version 1.5.1; celescope1.5.1b0

options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(patchwork)})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--mode", help = "mode")
argv <- add_argument(argv, "--anno_mode", help = "cloud anno, Default: F")
argv <- add_argument(argv, "--fj_path", help = "fj path, split by ,")
argv <- add_argument(argv, "--fj_spname", help = "fj spname, split by ,")
argv <- add_argument(argv, "--rna_spname", help = "rna spname, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: EBV_stat")
argv <- parse_args(argv)

rds <- argv$rds
mode <- argv$mode
anno_mode <- ifelse(is.na(argv$anno_mode), "anno", argv$anno_mode)
fj_path <- argv$fj_path
fj_spname <- unlist(strsplit(argv$fj_spname, split = ","))
fj_mtx <- paste0(fj_path, "/", fj_spname, "/03.assemble/count.txt")
rna_spname <- unlist(strsplit(argv$rna_spname, split = ","))
print(fj_mtx)

outdir <- argv$outdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

# readRDS ----
data_seurat <- readRDS(rds)
if(anno_mode == "clound"){
    data_seurat$sample <- data_seurat$`Sample ID`
    data_seurat$cluster <- data_seurat$annot_full
}
if(anno_mode == "rd_sr"){
    data_seurat$sample <- data_seurat$orig.ident
    data_seurat$cluster <- data_seurat$celltype
}

print(table(data_seurat$sample))

# cluster
options(repr.plot.height = 5, repr.plot.width = 7)
cluster_plot_out = lapply(unique(data_seurat$sample), function(x){
    DimPlot(subset(data_seurat, subset = sample == x), group.by = "cluster", label = T)
    ggsave(str_glue("{outdir}/cluster_{x}.pdf"), height = 5, width = 7)
    ggsave(str_glue("{outdir}/cluster_{x}.png"), height = 5, width = 7)
})

# sample_cluster
#options(repr.plot.height = 5, repr.plot.width = 5)
#data_seurat@meta.data %>%
#    ggplot(aes(x=sample, fill=cluster))+
#        geom_bar(color="black",position="fill",width = 0.5)+
#        theme_bw()+
#        theme(legend.title = element_blank(),axis.title = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))

print(table(data_seurat$cluster, data_seurat$sample))

# sample cluster table --
df_sample_cluster <- as.data.frame(t(as.data.frame(table(data_seurat$sample, data_seurat$cluster)) %>% spread(Var2, Freq)))
row.names(df_sample_cluster)[1] <- "cluster"

write.table(df_sample_cluster, str_glue("{outdir}/sample_cluster.xls"), quote=F, sep="\t", col.names=F)


# CR stat ----
function_read_fj_mtx <- function(df, spname){
    df <- read.table(df,sep="\t", header = T, row.names = 1)
    if(anno_mode!="rd_sr"){
        row.names(df) <- paste0(spname, "_", row.names(df))
    }
    df <- df[intersect(colnames(data_seurat),row.names(df)),]
    #print(dim(df))
    return(df)
}

df_input <- data.frame(fj_mtx= fj_mtx,
                       spname = rna_spname)

df_fj <- do.call(rbind, lapply(1:nrow(df_input), function(x){ function_read_fj_mtx(df_input[x,1], df_input[x,2]) }))

df_fj <- left_join(data.frame(bc = colnames(data_seurat)), 
                   rownames_to_column(df_fj, var="bc"), by = c("bc"="bc"))  # match zl
row.names(df_fj) <- df_fj$bc
df_fj[ is.na(df_fj)] <- 0

if(!identical(row.names(df_fj), colnames(data_seurat))){
    print("fail")
    quit()
}
print(identical(row.names(df_fj), colnames(data_seurat)))

if(mode == "TCR"){
    data_seurat$TCR_status <- factor(ifelse(df_fj$mark == "CB", "TCR+", "TCR-"), levels = c("TCR+", "TCR-"))
    data_seurat$cal_positive <- ifelse(df_fj$mark == "CB", 1, 0)
}

if(mode == "BCR"){
    data_seurat$BCR_status <- factor(ifelse(df_fj$mark == "CB", "BCR+", "BCR-"), levels = c("BCR+", "BCR-"))
    data_seurat$cal_positive <- ifelse(df_fj$mark == "CB", 1, 0)
}


function_stat <- function(spname, outdir, mode){
    data_subset <- subset(data_seurat, subset = sample == spname)
    spname <- gsub("RNA", mode, spname)   # trans zl to T(B)CR
    if(length(unique(data_subset@meta.data[, paste0(mode,"_status")]))==1){
        # just negativate
        DimPlot(data_subset, group.by = paste0(mode,"_status"), cols = c("lightgrey"))
        ggsave(str_glue("{outdir}/{spname}_{mode}_status.pdf"), height = 5, width = 7)
        ggsave(str_glue("{outdir}/{spname}_{mode}_status.png"), height = 5, width = 7)       
    }else{
        DimPlot(data_subset, group.by = paste0(mode,"_status"), cols = c("red","lightgrey"))
        ggsave(str_glue("{outdir}/{spname}_{mode}_status.pdf"), height = 5, width = 7)
        ggsave(str_glue("{outdir}/{spname}_{mode}_status.png"), height = 5, width = 7)    
    }

    stat_df <- data_subset@meta.data %>%
        group_by(cluster) %>%
        summarise(tag = mode,
                  positive_cells = sum(cal_positive >0),
                  total_cells = length(cal_positive)) %>%
        mutate(detect_percent = paste0(round(100*positive_cells/total_cells, 2), " %"))
    write.table(stat_df, str_glue("{outdir}/{spname}_{mode}_stat.tsv"), sep="\t", row.names = F, quote = F)
}

out = lapply(rna_spname, function_stat, outdir = outdir, mode = mode)