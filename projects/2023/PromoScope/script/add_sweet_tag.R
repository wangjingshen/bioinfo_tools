suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds of Main")
argv <- add_argument(argv,"--rds_fix", help="fix rds from h5ad")
argv <- add_argument(argv,"--sweet_cols", help="sweet_cols, v1 or v2, Default:v1")
argv <- add_argument(argv,"--cluster_cols", help="cluster_cols")
argv <- add_argument(argv,"--group_cols", help="group_cols")
argv <- add_argument(argv,"--tag_file", help="tag file, split by ,")
argv <- add_argument(argv,"--sname", help="sample name, split by ,")
argv <- add_argument(argv,"--subcluster_rds", help="Optional, rds of subcluster")
argv <- add_argument(argv,"--filter_cluster", help="Optional, filter cluster")
argv <- add_argument(argv,"--fill_group", help="Optional, fill group, default: F")
argv <- add_argument(argv,"--fill_sample", help="Optional, fill sample, default: F")
argv <- add_argument(argv,"--split_group", help="Optional, split group, default: F")
argv <- add_argument(argv,"--split_sample", help="Optional, split sample, default: F")
argv <- add_argument(argv,"--add_p", help="Optional, add_p, default: F")
argv <- add_argument(argv,"--outdir", help="outdir, such as: sweet_tag")
#argv <- add_argument(argv,"--name", help="name, project id + cluster name, such as: P23051601_Main")
argv <- add_argument(argv,"--saveRDS", help="saveRDS, default: T")
argv <- parse_args(argv)

rds <- argv$rds
rds_fix <- ifelse(is.na(argv$rds_fix), "F", argv$rds_fix)
sweet_cols <- ifelse(is.na(argv$sweet_cols), "v1", argv$sweet_cols)

if(!is.na(argv$cluster_cols)){
    cluster_cols <- unlist(str_split(argv$cluster_cols, ','))
}else{
    cluster_cols <- color_protocol
}
#print(cluster_cols)
if(!is.na(argv$group_cols)){
    group_cols <- unlist(str_split(argv$group_cols, ','))
}else{
    group_cols <- color_protocol
}
tag_file <- unlist(str_split(argv$tag_file, ','))
print(tag_file)
sname <- unlist(str_split(argv$sname, ','))
if(!is.na(argv$subcluster_rds)){
    subcluster_rds <- argv$subcluster_rds
}

fill_group <- ifelse(is.na(argv$fill_group), "F", argv$fill_group)
fill_sample <- ifelse(is.na(argv$fill_sample), "F", argv$fill_sample)
split_group <- ifelse(is.na(argv$split_group), "F", argv$split_group)
split_sample <- ifelse(is.na(argv$split_sample), "F", argv$split_sample)
add_p <- ifelse(is.na(argv$add_p), "F", argv$add_p)

outdir <- argv$outdir
#name <- argv$name
#res_outdir <- paste0(outdir, "/result/", name)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}
#if(!dir.exists(res_outdir)){
#    dir.create(res_outdir, recursive = T)
#}

if(is.na(argv$saveRDS)){
    saveRDS <- "F"
}else{
    saveRDS <- argv$saveRDS
}


#### function ----
# add sweet tag --
function_add_sweet_tag <- function(data, tag_file, prefix){
    function_read_sweet_tag <- function(i){
        tag <- read.table(tag_file[i], sep="\t", header = T, row.names = 1)
        tag$barcode <- paste0(prefix[i], "_", row.names(tag))
        row.names(tag) <- tag$barcode
        return(tag[,c("barcode","sweet_tag")])
    }
    tag_raw <- do.call(rbind, lapply(1:length(tag_file), function_read_sweet_tag))
    tag <- tag_raw[ colnames(data),]  # match 
    if(identical(tag$barcode, colnames(data)) & length(intersect(tag$barcode, colnames(data))) == ncol(data) & sum(is.na(tag$sweet_tag)) == 0){
        print("Barcodes match.")
        data$sweet_tag <- tag$sweet_tag
        data$log_sweet_tag <- log1p(data$sweet_tag)
        return(data)
    }else{
        print(length(intersect(tag$barcode, colnames(data))))
        print(length(intersect(tag_raw$barcode, colnames(data))))
        print(sum(is.na(tag$sweet_tag)))
        print("Barcodes mismatch.")
        quit()
    }
}

# featureplot --
function_featureplot <- function(data, outdir){
    if(sweet_cols == "v1"){
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c("lightgrey","red")) + 
            ggtitle("sweet_tag")
    }else{
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c('#1E90FF','#FF0000')) + 
            ggtitle("sweet_tag")
    }

    ggsave(paste0(outdir, "/Featureplot_sweet_tag.pdf"), plot = p, height = 5, width = 6)
    ggsave(paste0(outdir, "/Featureplot_sweet_tag.png"), plot = p, height = 5, width = 6)
}

function_featureplot_split_group <- function(data, subset_group, outdir){
    data <- subset(data, subset = group == subset_group)
    print(paste0("cellnumber of ", subset_group, ": ",  ncol(data)))
    if(sweet_cols == "v1"){
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c("lightgrey","red")) + 
            ggtitle(paste0("sweet_tag: ", subset_group))
    }else{
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c('#1E90FF','#FF0000')) + 
            ggtitle(paste0("sweet_tag: ", subset_group))
    }

    ggsave(paste0(outdir, "/Featureplot_sweet_tag_group_", subset_group, ".pdf"), plot = p, height = 5, width = 6)
    ggsave(paste0(outdir, "/Featureplot_sweet_tag_group_", subset_group, ".png"), plot = p, height = 5, width = 6)
}

function_featureplot_split_sample <- function(data, subset_sample, outdir){
    data <- subset(data, subset = sample == subset_sample)
    print(paste0("cellnumber of ", subset_sample, ": ",  ncol(data)))
    if(sweet_cols == "v1"){
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c("lightgrey","red")) + 
            ggtitle(paste0("sweet_tag: ", subset_sample))
    }else{
        p <- FeaturePlot(data, features = "sweet_tag_CLR", cols = c('#1E90FF','#FF0000')) + 
            ggtitle(paste0("sweet_tag: ", subset_sample))
    }

    ggsave(paste0(outdir, "/Featureplot_sweet_tag_sample_", subset_sample, ".pdf"), plot = p, height = 5, width = 6)
    ggsave(paste0(outdir, "/Featureplot_sweet_tag_sample_", subset_sample, ".png"), plot = p, height = 5, width = 6)
}


# vlnplot --
function_vlnplot <- function(data, outdir){
    filter_var <- names(table(data$cluster))[table(data$cluster) == 1]
    if(length(filter_var) >0 ){
        subset_var <- setdiff(unique(data$cluster), filter_var)
        data <- subset(data, subset = cluster %in% subset_var)  # subset
    }

    p <- data@meta.data %>%
        ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = cluster))+
            geom_violin(scale = 'width')+
            #geom_boxplot(width = 0.1,outlier.shape = NA)+
            guides(fill = guide_legend(title = "group"))+
            scale_fill_manual(values = cluster_cols)+
            ylab("sweet_tag")+
            #ggtitle("sweet_tag")+
            #ggpubr::stat_compare_means(label = "p.format")+
            theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                    plot.title = element_blank(),
                                    axis.title.x = element_blank(),
                                    legend.position = 'none')

    width = 3 + 0.5*length(unique(data$cluster))
    #print(width)
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag.pdf"), plot = p, height = 5, width = width)
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag.png"), plot = p, height = 5, width = width)
}

# boxplot --
function_boxplot <- function(data, outdir){
    data@meta.data %>%
        ggplot(aes(x= cluster, y = sweet_tag_CLR, fill = cluster)) +
            geom_boxplot(outlier.shape = NA) +
            scale_fill_manual(values = cluster_cols) +
            ylab("sweet_tag")+
            #ggtitle("sweet_tag") +
            theme_classic() +
            theme(axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
                legend.position = 'none',
                plot.title = element_blank())
    width = 3 + 0.5*length(unique(data$cluster))
    #print(width)
    ggsave(paste0(outdir, "/Boxplot_sweet_tag.pdf"), height = 5, width = width)
    ggsave(paste0(outdir, "/Boxplot_sweet_tag.png"), height = 5, width = width)
}


# stat df --
function_stat_df <- function(data, by, outdir){
    if(by=="cluster"){
        mean_df <- aggregate(data@meta.data$sweet_tag_CLR, by = list(data$cluster), FUN=mean)
        colnames(mean_df) <- c("cluster", "sweet_value")
    }
    if(by=="cluster_group"){
        mean_df <- aggregate(data@meta.data$sweet_tag_CLR, by = list(data$cluster, data$group), FUN=mean)
        colnames(mean_df) <- c("cluster", "group", "sweet_value")
    }
    if(by=="cluster_sample"){
        mean_df <- aggregate(data@meta.data$sweet_tag_CLR, by = list(data$cluster, data$sample), FUN=mean)
        colnames(mean_df) <- c("cluster", "sample", "sweet_value")
    }
    write.table(mean_df, file= str_glue("{outdir}/sweet_tag_mean_{by}.xls"), sep='\t', quote=F, row.names=F)
}

# split plot --
function_split <- function(data, split_var){
    data_split <- SplitObject(data, split.by = split_var)
    lapply(names(data_split), function(x){
        split_res_outdir <- paste0(outdir, "/split/", split_var, "_", x)
        if(!dir.exists(split_res_outdir)){
            dir.create(split_res_outdir, recursive = T)
        }
        function_featureplot(data_split[[x]], split_res_outdir)
        function_vlnplot(data_split[[x]], split_res_outdir)
        function_boxplot(data_split[[x]], split_res_outdir)
        function_stat_df(data_split[[x]], split_var, split_res_outdir) 
    })
}  


function_violin_boxplot_fill_group <- function(data, add_p, outdir){
    data$cluster_group <- paste0(data$cluster, "_", data$group)
    filter_var <- names(table(data$cluster_group))[table(data$cluster_group) == 1]
    if(length(filter_var) >0 ){
        subset_var <- setdiff(unique(data$cluster_group), filter_var)
        print(paste0("filter ", str_c(filter_var, collapse = " ")))
        data <- subset(data, subset = cluster_group %in% subset_var)  # subset
    }   
    if(add_p == "F"){
        p1 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
                geom_violin(scale = 'width')+
                #geom_boxplot(width = 0.2,position = position_dodge(width = 0.9), outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank())
        p2 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
                #geom_violin(scale = 'width')+
                geom_boxplot(outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank())
    }else{
        p1 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
                geom_violin(scale = 'width')+
                #geom_boxplot(width = 0.2,position = position_dodge(width = 0.9), outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                #ggtitle("sweet_tag")+
                ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank()) 
        p2 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
                #geom_violin(scale = 'width')+
                geom_boxplot(outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                #ggtitle("sweet_tag")+
                ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank()) 
    }
    width = 3 + 1*length(unique(data$cluster))
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag_cluster_group.pdf"), plot = p1, height = 5, width = width)
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag_cluster_group.png"), plot = p1, height = 5, width = width)
    ggsave(paste0(outdir, "/Boxplot_sweet_tag_cluster_group.pdf"), plot = p2, height = 5, width = width)
    ggsave(paste0(outdir, "/Boxplot_sweet_tag_cluster_group.png"), plot = p2, height = 5, width = width)    
}

function_violin_boxplot_fill_sample <- function(data, add_p, outdir){
    data$cluster_sample <- paste0(data$cluster, "_", data$sample)
    filter_var <- names(table(data$cluster_sample))[table(data$cluster_sample) == 1]
    if(length(filter_var) >0 ){
        subset_var <- setdiff(unique(data$cluster_sample), filter_var)
        print(paste0("filter ", str_c(filter_var, collapse = " ")))
        data <- subset(data, subset = cluster_sample %in% subset_var)  # subset
    }   
    if(add_p == "F"){
        p1 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = sample))+
                geom_violin(scale = 'width')+
                #geom_boxplot(width = 0.2,position = position_dodge(width = 0.9), outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank())
        p2 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = sample))+
                #geom_violin(scale = 'width')+
                geom_boxplot(outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank())
    }else{
        p1 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = sample))+
                geom_violin(scale = 'width')+
                #geom_boxplot(width = 0.2,position = position_dodge(width = 0.9), outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                #ggtitle("sweet_tag")+
                ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank()) 
        p2 <- data@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = sample))+
                #geom_violin(scale = 'width')+
                geom_boxplot(outlier.shape = NA)+
                guides(fill = guide_legend(title = "group"))+
                scale_fill_manual(values = group_cols)+
                ylab("sweet_tag")+
                #ggtitle("sweet_tag")+
                ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test")+
                theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                        plot.title = element_blank(),
                                        axis.title.x = element_blank()) 
    }
    width = 5 + 1.2*length(unique(data$cluster))
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag_cluster_sample.pdf"), plot = p1, height = 5, width = width)
    ggsave(paste0(outdir, "/Vlnplot_sweet_tag_cluster_sample.png"), plot = p1, height = 5, width = width)
    ggsave(paste0(outdir, "/Boxplot_sweet_tag_cluster_sample.pdf"), plot = p2, height = 5, width = width)
    ggsave(paste0(outdir, "/Boxplot_sweet_tag_cluster_sample.png"), plot = p2, height = 5, width = width)    
}


#### analysis ----
data_seurat <- readRDS(argv$rds)
if(rds_fix == "T"){
    data_seurat@assays$RNA@key <- "rna_"  # if h5ad to rds
    data_seurat$group <- data_seurat$gname
}
data_seurat <- function_add_sweet_tag(data_seurat, tag_file, sname)
# Create ADT
data_seurat[["ADT"]] <- CreateAssayObject(counts = t(data.frame(row.names = names(data_seurat$sweet_tag),
                                                     sweet = as.numeric(data_seurat$sweet_tag))))
data_seurat_list <- SplitObject(data_seurat, split.by = "sample")
data_seurat_list <- lapply(data_seurat_list, NormalizeData, assay = "ADT", margin = 1, normalization.method = "CLR")

data_seurat_merge <- merge(data_seurat_list[[1]], data_seurat_list[-1])   # merge CLR
# check merge and raw
print(paste0("all barcodes are equal: ", identical(colnames(data_seurat), colnames(data_seurat_merge))))
# merge之后没有reduction，这里将merge后的CLR加到原对象的metadata中
data_seurat$sweet_tag_CLR <- as.numeric(data_seurat_merge@assays$ADT@data[1,])
# 根据中值设置高低组
data_seurat$sweet_tag_group <- ifelse(data_seurat$sweet_tag_CLR > median(data_seurat$sweet_tag_CLR), "sweet_tag_high", "sweet_tag_low")
data_seurat$sweet_tag_group <- factor(data_seurat$sweet_tag_group, levels = c("sweet_tag_high", "sweet_tag_low"))
print(table(data_seurat$sweet_tag_group))
print(paste0("max value in sweet_tag_low: ", max(data_seurat$sweet_tag_CLR[data_seurat$sweet_tag_group == "sweet_tag_low"])))
print(paste0("min value in sweet_tag_high: ", min(data_seurat$sweet_tag_CLR[data_seurat$sweet_tag_group == "sweet_tag_high"])))


if(is.na(argv$subcluster_rds)){
    function_featureplot(data_seurat, outdir)
    function_vlnplot(data_seurat, outdir)
    function_boxplot(data_seurat, outdir)
    function_stat_df(data_seurat, by = "cluster", outdir)
    
    if(fill_group %in% c("T", "True", "TRUE")){
        lapply(unique(data_seurat$group), function_featureplot_split_group, data = data_seurat, outdir = outdir)
        function_violin_boxplot_fill_group(data_seurat, add_p = add_p, outdir)
        function_stat_df(data_seurat, by = "cluster_group", outdir)
    }

    if(fill_sample %in% c("T", "True", "TRUE")){
        lapply(unique(data_seurat$sample), function_featureplot_split_sample, data = data_seurat, outdir = outdir)
        function_violin_boxplot_fill_sample(data_seurat, add_p = add_p, outdir)
        function_stat_df(data_seurat, by = "cluster_sample", outdir)
    }

    if(saveRDS %in% c("T","True","TRUE")){
        saveRDS(data_seurat, paste0(outdir, "/sweet_tag.rds"))
    }

    df_sweet_tag <- data.frame( barcode = colnames(data_seurat),
                                sweet_tag_CLR = data_seurat@meta.data$sweet_tag_CLR)
    write.table(df_sweet_tag, paste0(outdir, "/sweet_tag.tsv"), sep = "\t", row.names = F, quote = F)


    if(split_group %in% c("T", "True", "TRUE")){
        function_split(data = data_seurat, "group")
        print("----Split group done.----")
    }
    if(split_sample %in% c("T", "True", "TRUE")){
        function_split(data = data_seurat, "sample")
        print("----Split sample done.----")
    }
}


# subcluster add sweet tag --
if(!is.na(argv$subcluster_rds)){
    # 读取分群,画featureplot
    data_subcluster <- readRDS(subcluster_rds)
    data_subset <- subset(data_seurat, cells = colnames(data_subcluster))
    # check
    print(paste0("subcluster barcodes are equal: ",identical(colnames(data_subcluster), colnames(data_subset))))
    data_subcluster$sweet_tag <- data_subset$sweet_tag
    data_subcluster$log_sweet_tag <- data_subset$log_sweet_tag 
    data_subcluster$sweet_tag_CLR <- data_subset$sweet_tag_CLR
    data_subcluster$sweet_tag_group <- data_subset$sweet_tag_group
    if(!is.na(argv$filter_cluster)){
        data_subcluster <- subset(data_subcluster, subset = cluster %in% setdiff(unique(data_subcluster$cluster), unlist(str_split(argv$filter_cluster, ','))))
    }

    function_featureplot(data_subcluster, outdir)
    function_vlnplot(data_subcluster, outdir)
    function_boxplot(data_subcluster, outdir)
    function_stat_df(data_subcluster, "cluster", outdir)
    
    if(fill_group %in% c("T", "True", "TRUE")){
        lapply(unique(data_subcluster$group), function_featureplot_split_group, data = data_subcluster, outdir = outdir)
        function_violin_boxplot_fill_group(data_subcluster, add_p=add_p, outdir)
        function_stat_df(data_subcluster, "cluster_group", outdir)
    }

    if(fill_sample %in% c("T", "True", "TRUE")){
        lapply(unique(data_subcluster$sample), function_featureplot_split_sample, data = data_subcluster, outdir = outdir)
        function_violin_boxplot_fill_sample(data_subcluster, add_p=add_p, outdir)
        function_stat_df(data_subcluster, "cluster_sample", outdir)
    }

    if(saveRDS %in% c("T","True","TRUE")){
        saveRDS(data_subcluster, paste0(outdir, "/sweet_tag_subcluster.rds"))
    }

    df_sweet_tag <- data.frame( barcode = colnames(data_subcluster),
                                sweet_tag_CLR = data_subcluster@meta.data$sweet_tag_CLR)
    write.table(df_sweet_tag, paste0(outdir, "/sweet_tag_subcluster.tsv"), sep = "\t", row.names = F, quote = F)

    if(split_group %in% c("T", "True", "TRUE")){
        function_split(data = data_subcluster, "group")
        print("----Split group done.----")
    }
    if(split_sample %in% c("T", "True", "TRUE")){
        function_split(data = data_subcluster, "sample")
        print("----Split sample done.----")
    }
}

print("----------")
print("Done.")
