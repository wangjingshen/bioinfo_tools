suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds of Main")
argv <- add_argument(argv,"--sweet_cols", help="sweet_cols, v1 or v2, Default:v1")
argv <- add_argument(argv,"--cluster", help="cluster, split by ,")
argv <- add_argument(argv,"--fill", help="group or sample")
argv <- add_argument(argv,"--outdir", help="outdir, such as: sweet_tag")
argv <- add_argument(argv,"--plot_box", help="plot box, yes or no. Default:no")
#argv <- add_argument(argv,"--name", help="name, project id + cluster name, such as: P23051601_Main")
argv <- parse_args(argv)

rds <- argv$rds
sweet_cols <- ifelse(is.na(argv$sweet_cols), "v1", argv$sweet_cols)

if(is.na(argv$sweet_cols)){
    cluster=""
}else{
    cluster <- unlist(str_split(argv$cluster, ','))
}
fill<- argv$fill
outdir <- argv$outdir
plot_box <- ifelse(is.na(argv$plot_box), "no", argv$plot_box)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = T)
}

#### function ----

# featureplot --
function_featureplot <- function(data, subset_group = "", subset_sample = "", outdir, name){  #subset_cluster = "",
    title = ""
    #if(subset_cluster != ""){
    #    print("cellnumber of data:")
    #    print(table(data$cluster))
    #    data <- subset(data, subset = cluster == subset_cluster)
    #    print(paste0("cellnumber of ", subset_cluster, ": ",  ncol(data)))
    #    title <- paste0(": ", subset_cluster)
    #}
    if(subset_group != ""){
        data <- subset(data, subset = group == subset_group)
        print(paste0("cellnumber of ", subset_group, ": ",  ncol(data)))
        title <- paste0(": ", subset_group)
    }
    if(subset_sample != ""){
        data <- subset(data, subset = sample == subset_sample)
        print(paste0("cellnumber of ", subset_sample, ": ",  ncol(data)))
        title <- paste0(": ", subset_sample)
    }
    plot_title <- paste0("sweet_tag", title, "_", name)
    if(sweet_cols == "v1"){
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c("lightgrey","red")) + 
            ggtitle(plot_title)
    }else{
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c('#1E90FF','#FF0000')) + 
            ggtitle(plot_title)
    }

    ggsave(paste0(outdir, "/Featureplot_sweet_tag_", gsub(": ", "", title), "_", name, ".pdf"), plot = p, height = 5, width = 6)
    ggsave(paste0(outdir, "/Featureplot_sweet_tag_", gsub(": ", "", title), "_", name, ".png"), plot = p, height = 5, width = 6)
}

# vlnplot --
function_vlnplot_boxplot <- function(data, fill, outdir, name){  #subset_cluster = "", 
    #if(subset_cluster == ""){
    #    filter_var <- names(table(data$cluster))[table(data$cluster) == 1]
    #    if(length(filter_var) >0 ){
    #        subset_var <- setdiff(unique(data$cluster), filter_var)
    #        data <- subset(data, subset = cluster %in% subset_var)  # subset
    #    }
    if(fill == "group"){
        p1 <- data@meta.data %>%
            ggplot(aes(x = group, y = sweet_tag_CLR, fill = group))+
                geom_violin(scale = 'width')+
                geom_boxplot(width = 0.1,outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = color_protocol)+
                ylab("sweet_tag")+
                ggtitle(paste0("sweet_tag: ", name))+
                #ggpubr::stat_compare_means(label = "p.format")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_text(hjust = 0.5),
                                        axis.title.x = element_blank(),
                                        legend.position = 'none')
        p2 <- data@meta.data %>%
            ggplot(aes(x= group, y = sweet_tag_CLR, fill = group)) +
                geom_boxplot(outlier.shape = NA) +
                scale_fill_manual(values = color_protocol) +
                ylab("sweet_tag")+
                ggtitle(paste0("sweet_tag: ", name)) +
                theme_classic() + theme(axis.title.x = element_blank(),
                                        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
                                        legend.position = 'none',
                                        plot.title = element_text(hjust = 0.5))
        width = 3 + 1*length(unique(data$group))
        ggsave(paste0(outdir, "/Vlnplot_sweet_tag_group_", name, ".pdf"), plot = p1, height = 5, width = width)
        ggsave(paste0(outdir, "/Vlnplot_sweet_tag_group_", name, ".png"), plot = p1, height = 5, width = width)
        if(plot_box == "yes"){
            ggsave(paste0(outdir, "/Boxplot_sweet_tag_group_", name, ".pdf"), plot = p2, height = 5, width = width)
            ggsave(paste0(outdir, "/Boxplot_sweet_tag_group_", name, ".png"), plot = p2, height = 5, width = width)
        }
    }
    if(fill == "sample"){
        p1 <- data@meta.data %>%
            ggplot(aes(x = sample, y = sweet_tag_CLR, fill = sample))+
                geom_violin(scale = 'width')+
                geom_boxplot(width = 0.1,outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = color_protocol)+
                ylab("sweet_tag")+
                ggtitle(paste0("sweet_tag: ", name))+
                #ggpubr::stat_compare_means(label = "p.format")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_text(hjust = 0.5),
                                        axis.title.x = element_blank(),
                                        legend.position = 'none')
        p2 <- data@meta.data %>%
            ggplot(aes(x= sample, y = sweet_tag_CLR, fill = sample)) +
                geom_boxplot(outlier.shape = NA) +
                scale_fill_manual(values = color_protocol) +
                ylab("sweet_tag")+
                ggtitle(paste0("sweet_tag: ", name)) +
                theme_classic() + theme(axis.title.x = element_blank(),
                                        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
                                        legend.position = 'none',
                                        plot.title = element_text(hjust = 0.5))
        width = 3 + 1*length(unique(data$sample))
        #print(width)
        ggsave(paste0(outdir, "/Vlnplot_sweet_tag_sample_", name, ".pdf"), plot = p1, height = 5, width = width)
        ggsave(paste0(outdir, "/Vlnplot_sweet_tag_sample_", name, ".png"), plot = p1, height = 5, width = width)
        if(plot_box == "yes"){
            ggsave(paste0(outdir, "/Boxplot_sweet_tag_sample_", name, ".pdf"), plot = p2, height = 5, width = width)
            ggsave(paste0(outdir, "/Boxplot_sweet_tag_sample_", name, ".png"), plot = p2, height = 5, width = width)
        }
    }
}

# stat df --
function_stat_df <- function(data, mod, outdir){
    if(mod == "cluster"){
        mean_df_Main <- aggregate(data@meta.data$sweet_tag_CLR, by = list(data$cluster), FUN=mean)
        colnames(mean_df_Main) <- c("cluster", "sweet_value")
    }
    if(mod == "cluster"){
        mean_df_Main <- aggregate(data@meta.data$sweet_tag_CLR, by = list(data$cluster, data$group), FUN=mean)
        colnames(mean_df_Main) <- c("cluster", "group", "sweet_value")
        write.table(mean_df_Main, file= paste0(outdir, "/sweet_tag_mean_cluster_group.xls"), sep='\t', quote=F, row.names=F)
    }

    write.table(mean_df_Main, file= paste0(outdir, "/sweet_tag_mean.xls"), sep='\t', quote=F, row.names=F)
}


#### analysis ----
data_seurat <- readRDS(argv$rds)
print("all data")
if(fill == "group"){
    lapply(unique(data_seurat$group), function_featureplot, data = data_seurat, subset_sample = "", outdir = outdir, name = "all")
    function_vlnplot_boxplot(data = data_seurat, fill = "group", outdir = outdir, name = "all")
}
if(fill == "sample"){
    lapply(unique(data_seurat$sample), function_featureplot, data = data_seurat, subset_group = "", outdir = outdir, name ="all")
    function_vlnplot_boxplot(data = data_seurat, fill = "sample", outdir = outdir, name = "all")
}



for( i in cluster){
    if(i == ""){
        break
    }else{
        data_cluster <- subset(data_seurat, subset = cluster == i)
        print(paste0("cluster: ", i))
        if(fill == "group"){
            lapply(unique(data_cluster$group), function_featureplot, data = data_cluster, subset_sample = "", outdir = outdir, name =  i)
            function_vlnplot_boxplot(data = data_cluster, fill = "group", outdir = outdir, name = i)
        }
        if(fill == "sample"){
            lapply(unique(data_cluster$sample), function_featureplot, data = data_cluster, subset_group = "", outdir = outdir, name =  i)
            function_vlnplot_boxplot(data = data_cluster, fill = "sample", outdir = outdir, name = i)
        }
    }
}


print("----------")
print("Done.")
