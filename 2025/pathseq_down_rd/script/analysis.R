# ref: Tumor microbiome links cellular programs and immunity in pancreatic cancer

suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
    library(circlize)
    library(vegan)
    library(clusterProfiler)
    library(org.Hs.eg.db)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--df_genus", help = "df_genus")
argv <- add_argument(argv, "--rna_spname", help = "rna sample name, split by ,")
argv <- add_argument(argv, "--spname", help = "16S sample name, split by ,")
argv <- add_argument(argv, "--umi_threshold", help = "genus detect threshold, Default: 2")
argv <- add_argument(argv, "--outdir", help = "output dir, default: outdir")
argv <- parse_args(argv)

rds <- argv$rds
df_genus <- unlist(strsplit(argv$df_genus, split = ","))
rna_spname <- unlist(strsplit(argv$rna_spname, split = ","))
spname <- unlist(strsplit(argv$spname, split = ","))
umi_threshold <- ifelse(is.na(argv$umi_threshold), 2, as.numeric(argv$umi_threshold))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
print("df_genus: ")
print(df_genus)
print("rna sample: ")
print(rna_spname)
print("sc16S sample: ")
print(spname)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
    dir.create(str_glue("{outdir}/01.base/"), recursive = T)
    dir.create(str_glue("{outdir}/02.diff/"), recursive = T)
    dir.create(str_glue("{outdir}/03./"), recursive = T)
    dir.create(str_glue("{outdir}/04./"), recursive = T)

}

# read rds --
data_seurat <- readRDS(rds)
data_seurat$barcode <- colnames(data_seurat)
print(table(data_seurat$sample))

# read genus --
args_df <- data.frame(df_genus = df_genus,
                      rna_spname = rna_spname)

df_genus <- lapply(1:nrow(args_df), function(x){
    df <- t(read.table(args_df[x,1], sep="\t", row.names = 1, header = T))
    row.names(df) <- paste0(args_df[x,2], "_", row.names(df))   
    return(df)
})
all_genus <- Reduce(union, lapply(df_genus, colnames))

# fill all_genus for each df_genus
df_genus <- lapply(df_genus, function(x){
    df <- data.frame(genus = all_genus) %>%
        left_join(rownames_to_column(as.data.frame(t(x)), var = "genus") , by = c("genus" = "genus")) %>%
        t() %>%
        as.data.frame()
    colnames(df) <- df[1, ]
    df <- df[-1, ]
    df[is.na(df)] <- 0
    df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.numeric)  # character to numeric
    return(df)
})
df_genus <- do.call(rbind, df_genus)
df_genus <- df_genus[ colnames(data_seurat),]

if(!identical(colnames(data_seurat), row.names(df_genus))){
    print("bc not identical, please check the input.")
    print(dim(data_seurat))
    print(dim(df_genus))
    print(colnames(data_seurat)[1:3])
    print(row.names(df_genus)[1:3])
}else{
    print("bc identical.")
    print(dim(data_seurat))
    print(dim(df_genus))
    print(colnames(data_seurat)[1:3])
    print(row.names(df_genus)[1:3])
}

data_seurat$total_genus <- rowSums(df_genus)
data_seurat@meta.data <- cbind.data.frame(data_seurat@meta.data, df_genus)

## base plot ----
# totol genus umi: cluster group --
options(repr.plot.height = 5, repr.plot.width =7)
data_seurat@meta.data %>% 
    group_by(cluster, group) %>%
    summarise(total_genus_umi = sum(total_genus)) %>%
    ggplot(aes(x=cluster, y=total_genus_umi, fill=group)) +
        geom_bar(color="black", stat = "identity", position = position_dodge()) +
        theme_classic() + 
        scale_fill_manual(values = color_protocol) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
ggsave(str_glue("{outdir}/01.base/total_genus_umi_barplot.png"), height = 5, width = 7)

# totol genus umi group --
options(repr.plot.height = 5, repr.plot.width =4)
data_seurat@meta.data %>% 
    group_by(group) %>%
    summarise(total_genus_umi = sum(total_genus)) %>%
    ggplot(aes(x=group, y=total_genus_umi, fill=group)) +
        geom_bar(color="black", stat = "identity", position = position_dodge(), width = 0.5) +
        theme_classic() + 
        scale_fill_manual(values = color_protocol)
ggsave(str_glue("{outdir}/01.base/total_genus_umi_group_barplot.png"), height = 5, width = 4)

# DimPlot genus_detect
data_seurat$genus_detect <- factor(ifelse(data_seurat$total_genus > umi_threshold, "Genus_detected", "No_genus_detected"), 
                                   levels = c("No_genus_detected","Genus_detected"))

DimPlot(data_seurat, group.by = "genus_detect", cols = c("lightgrey", "red"), pt.size = 0.2, order = T)
ggsave(str_glue("{outdir}/01.base/genus_detect_featureplot.png"), height = 5, width = 7)
# split by group
DimPlot(data_seurat, group.by = "genus_detect", cols = c("lightgrey", "red"), pt.size = 0.2, order = T, split.by = "group")
ggsave(str_glue("{outdir}/01.base/genus_detect_group_featureplot.png"), height = 5, width = 13)



## 02 diff ----
genus <- colnames(data_seurat@meta.data)[(which(colnames(data_seurat@meta.data) == "total_genus") + 1) :(which(colnames(data_seurat@meta.data) == "genus_detect") - 1)]
print("analysis_genus: ")
print(genus)
df_genus <- cbind.data.frame(data_seurat@meta.data[,genus], 
                             data_seurat@meta.data[,c("cluster","group")])

genus_diff_stat <- do.call(rbind,lapply(genus, function(x){
    data_seurat$analysis_genus <- data_seurat@meta.data[,x]
    stat_df <- ggpubr::compare_means(analysis_genus~group,data_seurat@meta.data,method = "wilcox.test") %>% dplyr::select(-.y.)
    mean_df <- aggregate(data_seurat$analysis_genus, by = list(data_seurat$group), FUN=mean)
    for (i in 1:nrow(mean_df)){
        stat_df[,paste0("mean_", mean_df[i,1])] = mean_df[i,2]
    }
    print("------")
    print(unlist(apply(stat_df[,grep("mean_", colnames(stat_df))], 1, which.max)))
    print(1:nrow(mean_df))
    print(colnames(stat_df)[grep("mean_", colnames(stat_df))])
    stat_df$mean_max <- plyr::mapvalues(apply(stat_df[,grep("mean_", colnames(stat_df))], 1, which.max), from = 1:nrow(mean_df), to = colnames(stat_df)[grep("mean_", colnames(stat_df))])
    stat_df$mean_max <- gsub("mean_", "", stat_df$mean_max)
    stat_df$genus <- x
    return(stat_df)
}))
#genus_diff_stat$genus <- genus
genus_diff_stat <- genus_diff_stat[ genus_diff_stat$p < 0.05,]
write.table(genus_diff_stat, str_glue("{outdir}/02.diff/stat.tsv"), sep="\t", quote=F, row.names = F)
genus_diff_sig <- genus_diff_stat$genus

for (i in genus_diff_sig){
    p1 <- VlnPlot(data_seurat, features = i, group.by = "group", split.by = "group", pt.size = 0, cols = color_protocol) + 
        ggpubr::stat_compare_means(label = "p.signif",method = "wilcox.test")
    ggsave(str_glue("{outdir}/02.diff/vlnplot_group_{i}.png"),plot =  p1, width = 7, height = 5)
}



saveRDS(data_seurat, "data_seurat.rds")


#df_genus_plot <- gather(df_genus, key = "genus", value = 'umi', -cluster, -group)
#function_diff_barplot <- function(g){
#    genus_diff_stat_subset <- genus_diff_stat[ , ]
#    genus_filter <- genus_diff_sig[ genus_diff_sig$group1 == g & genus_diff_sig$mean_max == g, ]$genus   
#    p <- df_genus_plot  %>%
#        filter(genus %in% genus_filter) %>%
#        group_by(genus, group) %>%
#        summarise(mean_umi = mean(umi)) %>%  #     summarise(total_umi = sum(umi)) %>%
#            ggplot(aes(group, mean_umi, fill = genus))+
#                geom_bar(color="black", stat="identity", width = 0.5)+
#                theme_classic()+
#                theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),axis.title.x = element_blank()) + 
#                scale_fill_manual(values = color_protocol[-13])   
#    return(p)
#}

#for (i in unique(data_seurat$group)){
#    function_diff_barplot(i)
#    ggsave(str_glue("{outdir}/02.diff/diff_genus_barplot_{i}.png"), height = 5, width = 5)
#}

