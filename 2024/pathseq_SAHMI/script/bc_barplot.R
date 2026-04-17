# down analysis for /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16s_sgr/script/pipeline.py
# top5的并集以外的当作others

suppressMessages(suppressWarnings({
    library(Seurat)
    library(argparser)
    library(tidyverse)
}))

N <- 10
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds, split by ,")
argv <- add_argument(argv, "--subset", help = "subset, split by ,(spname) and ;(rds)")
argv <- add_argument(argv, "--anno", help = "anno")
argv <- add_argument(argv, "--mode", help = "mode, Default: detect ")
argv <- add_argument(argv, "--plot_mode", help = "plot mode, Default: all,  ")
argv <- add_argument(argv, "--df_genus", help = "genus_umi.tsv, file form total_bacteria, split by , ")
argv <- add_argument(argv, "--rna_spname", help = "rna_spname, split by ,")
argv <- add_argument(argv, "--spname", help = "spname, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir, Default: outdir")
argv <- parse_args(argv)

rds <- unlist(str_split(argv$rds, pattern = ",", simplify = F))
df_genus <- unlist(str_split(argv$df_genus, pattern = ",", simplify = F))
rna_spname <- unlist(str_split(argv$rna_spname, pattern = ",", simplify = F))
spname <- unlist(str_split(argv$spname, pattern = ",", simplify = F))
mode <- ifelse(is.na(argv$mode), "detect", argv$mode)
plot_mode <- ifelse(is.na(argv$plot_mode), "all", argv$plot_mode)
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}


if(length(rds)==1){
    data_seurat <- readRDS(rds)
}else{
    data_seurat <- lapply(rds, readRDS)
    data_seurat <- merge(data_seurat[[1]], y = data_seurat[-1])
    #print(dim(data_seurat))
}

if(!is.na(argv$subset)){
    subset_list <- unlist(str_split(argv$subset, pattern = ",", simplify = F))
    print(table(data_seurat$sample))
    data_seurat <- subset(data_seurat, subset = sample %in% subset_list)
    print("filter:")
    print(table(data_seurat$sample))
}

if(argv$anno == "singleR"){
    data_seurat$cluster = data_seurat$cluster_singleR
}


df <- data.frame(spname = spname,
                 df_genus = df_genus,
                 rna_spname = rna_spname)

df_list <- lapply(1:nrow(df),function(x){
    df_genus <- read.table(df[x,2], header = T, row.names = 1, sep="\t")
    return(df_genus)
})

# 所有样本的genus
all_genus <- unique(unlist(lapply(df_list, function(x){    # 
    colnames(x)
})))
print(paste("all genus: ",length(all_genus)))

plot_data <- do.call(rbind,lapply(df_list, function(x){
    #print(dim(x))
    #print(length(setdiff(all_genus, colnames(x))))
    x[, setdiff(all_genus, colnames(x))] = 0
    #print(dim(x))
    #return(x)
    return(x[,all_genus])
}))
#print(row.names(plot_data)[1:5])
#print(colnames(data_seurat)[1:5])

if(identical(row.names(plot_data), colnames(data_seurat)) == F){
    print("bc not equal.")
    print(row.names(plot_data)[1:5])
    print(colnames(data_seurat)[1:5])
    quit()
}else{
    print("bc equal.")
}

plot_data$cluster <- data_seurat$cluster
plot_data$spname <- plyr::mapvalues(data_seurat$sample, from = rna_spname, to = spname)
plot_data <- gather(plot_data, key = "genus", value = 'umi', -spname, -cluster)

if(mode == "detect"){
    out <- lapply(unique(plot_data$cluster), function(x){
        plot_data_s <- plot_data[ plot_data$cluster == x & plot_data$umi>0, ] 
        #print(plot_data_s[1:2,])
        width =9
        if(plot_mode == "top"){
            tmp <- plot_data_s %>% group_by(spname, genus) %>%
                   summarise(detect_bc = length(umi),
                             mean_umi = mean(umi)) %>%
                   top_n(N, mean_umi) %>% do (head(., n=N))
                   #slice_max(mean_umi, n = N)
            top_genus <- unique(tmp$genus)
            plot_data_s[ !plot_data_s$genus %in% top_genus, "genus"] <- "others"
            plot_data_s$genus <- factor(plot_data_s$genus, levels = c(top_genus, "others"))
            width = 7
        }

        plot_data_s <- plot_data_s %>% group_by(spname, genus) %>%
                summarise(mean_umi = mean(umi))
        
        write.table(plot_data_s, file = str_glue('{outdir}/{x}_stat_detect.tsv'), sep="\t", quote=F, row.names = F )

        #write.table(tmp, file = str_glue('{outdir}/{x}_stat_detect_top5.tsv'), sep="\t", quote=F, row.names = F )

        options(repr.plot.width = 9, repr.plot.height = 6)
        ggplot(plot_data_s, aes(spname, mean_umi, fill=genus))+
            geom_bar(color="black", stat="identity", width = 0.2 )+
            scale_fill_manual(values = color_protocol)+
            ggtitle(x)+
            ylab("Average UMI per cell")+
            theme_classic()+
            theme(legend.position = 'right',
                  legend.title = element_blank(), 
                  plot.title = element_text(hjust = 0.5),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) 
        ggsave(str_glue('{outdir}/{x}_barplot_detect.png'), width = width, height = 6)
        ggsave(str_glue('{outdir}/{x}_barplot_detect.pdf'), width = width, height = 6)            
    })
}