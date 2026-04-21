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
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981",
                    "#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F",
                    "#000101","OrangeRed","SlateBlue","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink","Red",
                    "#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80",
                    "#483D8B","#66CD00","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23",
                    "#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",
                    "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726",
                    "#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD",
                    "#548B54","#BBFFFF","#B4EEB4","#00C5CD","#7FFFD4","#8EE5EE","#68838B","#B9D3EE","#9ACD32","#00688B","#FFEC8B",
                    "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#CD661D","#CDC5BF","#FF8C69",
                    "#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32",
                    "#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#191970","#E8E8E8")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--df_genus", help = "df_genus")
argv <- add_argument(argv, "--rna_spname", help = "rna sample name, split by ,")
argv <- add_argument(argv, "--spname", help = "sample name, split by ,")
argv <- add_argument(argv, "--subset", help = "subset sample, split by ,")
argv <- add_argument(argv, "--barplot_topn", help = "barplot_topn, default:5")
argv <- add_argument(argv, "--barplot_width", help = "barplot_width, default:7")
argv <- add_argument(argv, "--plot_total_barplot", help = "plot_total_barplot, default:F")
argv <- add_argument(argv, "--outdir", help = "output dir, default: outdir")
argv <- add_argument(argv, "--saveRDS", help = "saveRDS, default: F")
argv <- parse_args(argv)

rds <- argv$rds
df_genus <- unlist(strsplit(argv$df_genus, split = ","))
rna_spname <- unlist(strsplit(argv$rna_spname, split = ","))
spname <- unlist(strsplit(argv$spname, split = ","))
barplot_topn <- ifelse(is.na(argv$barplot_topn), 5, as.numeric(argv$barplot_topn))
barplot_width <- ifelse(is.na(argv$barplot_width), 7, as.numeric(argv$barplot_width))
plot_total_barplot <- ifelse(is.na(argv$plot_total_barplot), "F", as.numeric(argv$plot_total_barplot))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
saveRDS <- ifelse(is.na(argv$saveRDS), "F", argv$saveRDS)
print("df_genus: ")
print(df_genus)
print("rna sample: ")
print(rna_spname)
print("sc16S sample: ")
print(spname)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# read rds
data_seurat <- readRDS(rds)
if(!is.na(argv$subset)){
    data_seurat <- subset(data_seurat, subset = sample %in% unlist(strsplit(argv$subset, split = ",")))
}

data_seurat$barcode <- colnames(data_seurat)
print(table(data_seurat$sample))
#print(table(data_seurat$sample, data_seurat$cluster))

if(!dir.exists(str_glue("{outdir}/cluster/"))){
    dir.create(str_glue("{outdir}/cluster/"), recursive = T)
}
p <- DimPlot(data_seurat, group.by = "cluster", cols = color_protocol)
ggsave(str_glue("{outdir}/cluster/cluster.png"), plot = p, width = 7, height = 5)
# sample cluster table 
write.table(t(as.data.frame(table(data_seurat$sample, data_seurat$cluster)) %>% spread(Var2, Freq) ),
    str_glue("{outdir}/cluster/sample_cluster.tsv"), quote=F, sep="\t", col.names=F)

# read genus --
df_genus <- lapply(1:length(df_genus), function(x){
    df <- read.table(df_genus[x], sep="\t", row.names = 1, header = T)
    row.names(df) <- gsub(spname[x], rna_spname[x], row.names(df))   # trans df_genus spname to rna_spname
    #row.names(df) <- paste0(rna_spname[x], "_", row.names(df))  # for old version df genus, not use
    return(df)
})
all_genus <- Reduce(union, lapply(df_genus, colnames))

# merge all_genus
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
print(row.names(df_genus)[1:5])
print(colnames(data_seurat)[1:5])
df_genus <- df_genus[ match(colnames(data_seurat), row.names(df_genus)),]   # for seurat filter cells

# check barcode between seurat and df_genus
if( !identical(colnames(data_seurat), row.names(df_genus)) ){
    print("barcode not equal.")
    print(colnames(data_seurat)[1:3])
    print(row.names(df_genus)[1:3])
    print(dim(data_seurat))
    print(dim(df_genus))
    quit()
}else{
    print("barcode equal.")
}

# plot total barplot
if(plot_total_barplot == "T"){
    plot_data <- df_genus
    plot_data$cluster <- data_seurat$cluster
    plot_data$spname <- plyr::mapvalues(data_seurat$sample, from = rna_spname, to = spname) # trans zl sample to fj sample
    plot_data <- gather(plot_data, key = "genus", value = 'umi', -spname, -cluster)

    plot_data <- plot_data[ plot_data$umi>0, ]
    tmp <- plot_data_subset %>% group_by(spname, genus) %>%
        summarise(detect_bc = length(umi),
                  mean_umi = mean(umi)) %>%
        top_n(barplot_topn, mean_umi) %>% 
        do(head(., n=barplot_topn))   # some top_n values equal, selcet 'top'(head) barplot_topn
    top_spname_genus <- unique(paste0(tmp$spname, "_", tmp$genus))
    plot_data_subset$spname_genus <- paste0(plot_data_subset$spname, "_", plot_data_subset$genus)
    plot_data_subset[ !plot_data_subset$spname_genus %in% top_spname_genus, "genus"] <- "others"

    plot_data_subset <- plot_data_subset %>% 
        group_by(spname, genus) %>%
        summarise(mean_umi = mean(umi))
    
    # genus order umi, set levels
    plot_data_subset$genus <- factor(plot_data_subset$genus, 
        levels = c(setdiff(unique(plot_data_subset$genus[order(plot_data_subset$mean_umi, decreasing=T)]), "others"), "others"))

    # df
    write.table(plot_data_subset, file = str_glue('{outdir}/cluster_genus_mean_umi/{x}_genus_mean_umi.tsv'), sep="\t", quote=F, row.names = F )
    # plot
    ggplot(plot_data_subset, aes(spname, mean_umi, fill=genus))+
        geom_bar(color="black", stat="identity", width = 0.2)+
        ggtitle(x)+
        scale_fill_manual(values = color_protocol)+
        ylab("Average UMI per cell")+
        theme_classic()+
        theme(legend.title = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(),
              axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))
    ggsave(str_glue('{outdir}/cluster_genus_mean_umi/{x}_genus_mean_umi_barplot.png'), width = barplot_width, height = 6)
}


# add total_genus
data_seurat$total_genus <- rowSums(df_genus)
#
data_seurat@meta.data <- cbind.data.frame(data_seurat@meta.data, df_genus)


for (i in 1:length(rna_spname)){
    spname_subset <- spname[i]
    rna_spname_subset <- rna_spname[i]
    data_seurat_subset <- subset(data_seurat, subset = sample == rna_spname_subset)

    ## cluster
    DimPlot(data_seurat_subset, group.by = "cluster", cols = color_protocol)
    ggsave(str_glue("{outdir}/cluster/{spname_subset}_cluster.png"), width = 7, height = 5)

    ## total genus
    if(!dir.exists(str_glue("{outdir}/total_genus/"))){
        dir.create(str_glue("{outdir}/total_genus/"), recursive = T)
    }

    # plot total
    FeaturePlot(data_seurat_subset, 'total_genus', order = T, pt.size = 1)
    ggsave(str_glue('{outdir}/total_genus/{spname_subset}_total_genus_Featureplot.png'), width = 7, height = 5)

    # stat total
    stat_bc <- data_seurat_subset@meta.data %>%
        group_by(sample, cluster) %>%
        summarise(detect_bc = sum(total_genus > 0),
                  total_bc = length(total_genus),
                  mean_umi_in_detect_bc = sum(total_genus)/sum(total_genus>0) ) %>%
        mutate(detect_percent = paste0(round(100*detect_bc/total_bc, 2), " %"))
    write.table(stat_bc, str_glue("{outdir}/total_genus/{spname_subset}_total_genus_detect.tsv"),sep="\t", row.names = F, quote = F)

    ## top10 genus
    if(!dir.exists(str_glue("{outdir}/top10_genus/{spname_subset}/"))){
        dir.create(str_glue("{outdir}/top10_genus/{spname_subset}/"), recursive = T)
    }

    print("1")
    print(rna_spname_subset)
    print(row.names(df_genus)[1:5])
    print(grep(rna_spname_subset, row.names(df_genus))[1:2])
    df_genus_subset <- df_genus[ grep(paste0(rna_spname_subset, "_"), row.names(df_genus)),]  # get subset df_genus
    print(dim(df_genus_subset))
    print(dim(data_seurat_subset@meta.data))
    data_seurat_subset@meta.data <- cbind.data.frame(data_seurat_subset@meta.data, df_genus_subset)
    print("2")
    top10_genus <- names(sort(colSums(df_genus_subset > 0), decreasing = T)[1:10])
    print("top genus: ")
    print(top10_genus)
    top10_genus <- top10_genus[ !is.na(top10_genus)]  # filter NA, for top_n < 10
    print("top genus filter NA: ")
    print(top10_genus)

    # plot top10
    options(repr.plot.width=16, repr.plot.height = 10)
    FeaturePlot(data_seurat_subset, top10_genus, pt.size=1, order = T)
    ggsave(str_glue("{outdir}/top10_genus/{spname_subset}/{spname_subset}_top10_genus_detect_Featureplot.png"), width=16, height = 10)

    # split plot
    seq= 1
    for (genus_plot in top10_genus){
        FeaturePlot(data_seurat_subset, features = genus_plot, pt.size=1, order = T)
        ggsave(str_glue("{outdir}/top10_genus/{spname_subset}/{spname_subset}_top{seq}_{genus_plot}_Featureplot.png"), width=7, height = 5)
    
        data_seurat_subset$analysis_col <- data_seurat_subset@meta.data[,genus_plot]
        df_stat <- data_seurat_subset@meta.data %>% 
            group_by(cluster) %>%
            summarise(genus = genus_plot,
                      detect_bc = sum(analysis_col >0),
                      total_bc = length(analysis_col),
                      mean_umi_in_detect_bc = sum(analysis_col)/sum(analysis_col>0) ) %>%
            mutate(detect_percent = paste0(round(100*detect_bc/total_bc,2), " %"))  
        df_stat[is.na(df_stat)] <- 0
        write.table(df_stat, str_glue("{outdir}/top10_genus/{spname_subset}/{spname_subset}_top{seq}_{genus_plot}_detect.tsv"),sep="\t",row.names=F,quote = F)
        seq = seq + 1
    }

    if(saveRDS == "T"){
        if(!dir.exists(str_glue("{outdir}/rds/"))){
            dir.create(str_glue("{outdir}/rds/"), recursive = T)
        }
        saveRDS(data_seurat_subset, str_glue("{outdir}/rds/data_seurat_16S_{spname_subset}.rds"))
    }

}

## cluster genus mean UMI barplot 
if(!dir.exists(str_glue("{outdir}/cluster_genus_mean_umi/"))){
    dir.create(str_glue("{outdir}/cluster_genus_mean_umi/"), recursive = T)
}

plot_data <- df_genus
plot_data$cluster <- data_seurat$cluster
plot_data$spname <- plyr::mapvalues(data_seurat$sample, from = rna_spname, to = spname) # trans zl sample to fj sample
plot_data <- gather(plot_data, key = "genus", value = 'umi', -spname, -cluster)

out <- lapply(unique(data_seurat$cluster), function(x){
    plot_data_subset <- plot_data[ plot_data$cluster == x & plot_data$umi>0, ]
    #saveRDS(plot_data_subset, str_glue("{outdir}/{x}_plot_data.rds"))   # check
    tmp <- plot_data_subset %>% 
        group_by(spname, genus) %>%
        summarise(detect_bc = length(umi),
                  mean_umi = mean(umi)) %>%
        arrange(spname, desc(mean_umi)) %>%
        do(head(., n=barplot_topn))   
        #slice_max(mean_umi, n = N)
    #saveRDS(tmp, str_glue("{outdir}/{x}.rds"))  # check
    top_spname_genus <- unique(paste0(tmp$spname, "_", tmp$genus))
    plot_data_subset$spname_genus <- paste0(plot_data_subset$spname, "_", plot_data_subset$genus)
    plot_data_subset[ !plot_data_subset$spname_genus %in% top_spname_genus, "genus"] <- "others"
    #top_genus <- unique(tmp$genus)     # raw rule, sample may contain genus from others
    #plot_data_subset[ !plot_data_subset$genus %in% top_genus, "genus"] <- "others"
    #plot_data_subset$genus <- factor(plot_data_subset$genus, levels = c(top_genus, "others"))

    plot_data_subset <- plot_data_subset %>% 
        group_by(spname, genus) %>%
        summarise(mean_umi = mean(umi))
    
    # genus order umi, set levels
    plot_data_subset$genus <- factor(plot_data_subset$genus, 
        levels = c(setdiff(unique(plot_data_subset$genus[order(plot_data_subset$mean_umi, decreasing=T)]), "others"), "others"))

    # df
    write.table(plot_data_subset, file = str_glue('{outdir}/cluster_genus_mean_umi/{x}_genus_mean_umi.tsv'), sep="\t", quote=F, row.names = F )
    # plot
    ggplot(plot_data_subset, aes(spname, mean_umi, fill = genus))+
        geom_bar(color="black", stat="identity", width = 0.2)+
        ggtitle(x)+
        scale_fill_manual(values = color_protocol)+
        ylab("Average UMI per cell")+
        theme_classic()+
        theme(legend.title = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(),
              axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))

    ggsave(str_glue('{outdir}/cluster_genus_mean_umi/{x}_genus_mean_umi_barplot.png'), width = barplot_width, height = 6)
})

if(saveRDS == "T"){
    saveRDS(data_seurat, str_glue("{outdir}/rds/data_seurat_16S.rds"))
}


print('Done')