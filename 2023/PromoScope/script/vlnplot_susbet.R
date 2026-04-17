suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="sweet rds")
argv <- add_argument(argv,"--subset_group", help="subset group")
argv <- add_argument(argv,"--filter_cluster", help="Optional, filter cluster, default:F")
argv <- add_argument(argv,"--add_p", help="Optional, add_p, default: T")
argv <- add_argument(argv,"--outdir", help="outdir, such as: sweet_tag")
argv <- parse_args(argv)

rds <- argv$rds
subset_group <- unlist(str_split(argv$subset_group, ','))
filter_cluster <- ifelse(is.na(argv$filter_cluster), "F", unlist(str_split(argv$filter_cluster, ',')))
add_p <- ifelse(is.na(argv$add_p), "T", argv$add_p)
outdir <- argv$outdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

#### function ----
function_violin_boxplot_fill_group <- function(data, add_p, outdir, name){
    data$cluster_group <- paste0(data$cluster, "_", data$group)
    filter_var <- names(table(data$cluster_group))[table(data$cluster_group) == 1]
    if(length(filter_var) >0 ){
        subset_var <- setdiff(unique(data$cluster_group), filter_var)
        print(paste0("filter ", str_c(filter_var, collapse = " ")))
        data <- subset(data, subset = cluster_group %in% subset_var)  # subset
    }   
    if(add_p == "F"){
        p <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
                geom_violin(scale = 'width')+
                geom_boxplot(width = 0.2,position = position_dodge(width = 0.9), outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = color_protocol)+
                ylab("sweet_tag")+
                #ggtitle("sweet_tag")+
                #ggpubr::stat_compare_means(label = "p.format")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank())
    }else{
        p <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
                geom_violin(scale = 'width')+
                geom_boxplot(width = 0.2,position = position_dodge(width = 0.9), outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = color_protocol)+
                ylab("sweet_tag")+
                #ggtitle("sweet_tag")+
                ggpubr::stat_compare_means(label = "p.format", method = "wilcox.test")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank()) 
    }
    width = 3 + 1*length(unique(data$cluster))
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag_", name, ".pdf"), height = 5, width = width)
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag_", name, ".png"), height = 5, width = width)
    
}

# stat df --
function_stat_df_cluster_group <- function(data, outdir, name){
    mean_df_Main <- aggregate(data@meta.data$sweet_tag_CLR, by = list(data$cluster, data$group), FUN=mean)
    colnames(mean_df_Main) <- c("cluster", "group", "sweet_value")
    write.table(mean_df_Main, file= paste0(outdir, "/sweet_tag_mean_", name, ".xls"), sep='\t', quote=F, row.names=F)
}

#### analysis ----
data_seurat <- readRDS(argv$rds)
for (i in subset_group){
    i <- unlist(str_split(i, ':'))
    data_subset <- subset(data_seurat, subset = group %in% i)
    if(!is.na(argv$filter_cluster)){
        data_subset <- subset(data_subset, subset = cluster %in% setdiff(unique(data_subset$cluster), filter_cluster))
    }
    function_violin_boxplot_fill_group(data_subset, add_p = add_p, outdir, paste(i, collapse = "_vs_"))
    function_stat_df_cluster_group(data_subset, outdir, paste(i, collapse = "_vs_"))
}


print("----------")
print("Done.")
