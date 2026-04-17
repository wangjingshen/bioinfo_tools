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
argv <- add_argument(argv, "--mod", help = "same or no. Deafault: no")
argv <- add_argument(argv, "--anno_file", help = "anno file, for un-anno rds")
argv <- add_argument(argv, "--version", help = "version. v1 or v2. v2 for lims. Default: v2")
argv <- add_argument(argv, "--cele_version", help = "cele_version. Default: 1.5.1")
argv <- add_argument(argv, "--fj_path", help = "fj path, split by ,")
argv <- add_argument(argv, "--fj_name", help = "fj name, split by ,")
argv <- add_argument(argv, "--prefix", help = "prefix of cell name of zl, split by ,")
argv <- add_argument(argv, "--subset", help = "subset sample, split by ,")
argv <- add_argument(argv, "--sort", help = "plot cells in order of expression. Default:F")
#argv <- add_argument(argv, "--clean_name", help = "clean name, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: EBV_stat")
argv <- parse_args(argv)

# default value -- 
if(is.na(argv$version)){
    argv$version ="v2"
}
if(is.na(argv$mod)){
    argv$mod ="no"
}
if(is.na(argv$outdir)){
    argv$outdir ="EBV_stat/"
}
if(is.na(argv$cele_version)){
    argv$cele_version ="1.12"
}
argv$order <- ifelse(is.na(argv$sort), "F", argv$sort)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)   #  xxx_EBV_stat/
}

# version = argv$version # version is R function, not set, use argv$version
print(paste0("EBV stat version: ", argv$version))  # print version 

mod = argv$mod
outdir = argv$outdir
# fj_path 
if(argv$mod == "same"){
    fj_name <- unlist(strsplit(argv$fj_name, split = ","))
    fj_path <- paste0(argv$fj_path, "/", fj_name)
}
if(argv$mod == "no"){
    fj_path <- unlist(strsplit(argv$fj_path, split = ","))
    fj_name <- unlist(strsplit(argv$fj_name, split = ","))
}
print(fj_path)
print(fj_name)

if(!is.na(argv$prefix)){
    prefix = unlist(strsplit(argv$prefix, split = ","))
}
#clean_name = unlist(strsplit(argv$clean_name, split = ","))


# read rds -- 
data_seurat <- readRDS(argv$rds)
# anno for single un-anno rds, analysis as v1  ######## !!!!!!!!! be careful   !!!!!!!!!
if(!is.na(argv$anno_file)){
    anno_file <- read.table(argv$anno_file, header=T,sep="\t")
    anno_file <- anno_file[!duplicated(anno_file$cluster),]    # if duplicated, use first anno
    anno_file$cell_type[is.na(anno_file$cell_type)]="NA"     # to character, for plot label "NA"; if NA, not labeled
    if(length(intersect(unique(data_seurat$seurat_clusters), unique(anno_file$cluster)))){
        print("anno file cluster differernt seurat cluster")
        quit()
    }else{
        data_seurat@meta.data%>%
            left_join(anno_file,by = c("seurat_clusters" = "cluster")) -> anno_df
        data_seurat$cell_type <- anno_df$cell_type
        prefix = unique(data_seurat$orig.ident)   # fix auto assign 
    }
}

if(argv$version == "v1"){
    data_seurat$sample = data_seurat$orig.ident   # add sample col
    data_seurat$cluster = data_seurat$cell_type   # add cluster col
}

if(argv$version == "singleR"){
    data_seurat$sample = data_seurat$orig.ident   # add sample col
    data_seurat$cluster = data_seurat$celltype  # add cluster col
}

if(!is.na(argv$subset)){
    subset_sample <- unlist(strsplit(argv$subset, split = ","))
    data_seurat <- subset(data_seurat, subset = sample %in% subset_sample)   
}

function_read_fj_csv <- function(i){
    path = fj_path[i]
    if(argv$cele_version=="1.5.1"){
        fj <- read.table(paste0(path,"/06.analysis_capture_virus/",fj_name[i],"_virus_tsne.tsv"),sep="\t",header = T,row.names = 1)
        fj$UMI[ is.na(fj$UMI)] =0
        if(!is.na(argv$prefix)){
            bc_prefix = prefix[i]
            fj$barcode <- paste0(bc_prefix,"_",fj$barcode)
        }
        row.names(fj) <- fj$barcode
        return(fj[,c("barcode","UMI")])
    }else if(argv$cele_version=="1.1.8"){
        fj <- read.table(paste0(path,"/05.analysis_capture_virus/",fj_name[i],"_tsne.tsv"),sep="\t",header = T,row.names = 1)
        fj$UMI[ is.na(fj$UMI)] =0
        if(!is.na(argv$prefix)){
            bc_prefix = prefix[i]
            fj$barcode <- paste0(bc_prefix,"_",fj$barcode)
        }
        row.names(fj) <- fj$barcode
        return(fj[,c("barcode","UMI")])
    }else if(argv$cele_version=="1.7.2"){
        fj <- read.table(paste0(path,"/06.filter_virus/",fj_name[i],"_filtered_UMI_tsne.csv"),sep=",",header = T)
        if(!is.na(argv$prefix)){
            bc_prefix = prefix[i]
            fj$barcode <- paste0(bc_prefix,"_",fj$barcode)
        }
        row.names(fj) <- fj$barcode
        return(fj[,c("barcode","sum_UMI")])
    }else{
        fj <- read.table(paste0(path,"/07.analysis_virus/",fj_name[i],"_UMI_tsne.csv"),sep=",",header = T)
        if(!is.na(argv$prefix)){
            bc_prefix = prefix[i]
            fj$barcode <- paste0(bc_prefix,"_",fj$barcode)
        }
        row.names(fj) <- fj$barcode
        return(fj[,c("barcode","sum_UMI")])
    }

}
fj <- do.call(rbind,lapply(1:length(fj_name),function_read_fj_csv))
if(length(intersect(colnames(data_seurat),fj$barcode))!=ncol(data_seurat)){
    print("zl and fj not consist.")
    quit()
}else{
    fj <- fj[ colnames(data_seurat),]  # match 
    data_seurat$EBV_UMI <- fj[,2]    # support 1.13.0;  raw code: data_seurat$HBV_UMI <- fj$UMI   #data_seurat$HBV_UMI <- fj$sum_UMI
    data_seurat$EBV_status <- ifelse(data_seurat$EBV_UMI ==0, "EBV-", "EBV+")
    data_seurat$EBV_status <- factor(data_seurat$EBV_status,levels = c("EBV+","EBV-"))
}

# EBV stat --
function_EBV_plot_update <- function(i){
    if(!is.na(argv$prefix)){
        outfile_name <- name_df[i,1]
    }
    if(is.na(argv$prefix)){
        outfile_name <- name_df[i,2]
    }
    plot_data = data_seurat@reductions$umap@cell.embeddings %>% 
      as.data.frame() %>% 
      cbind(barcode = colnames(data_seurat))
    plot_data <- cbind(plot_data, data_seurat@meta.data)
    if(is.na(argv$prefix)){
        plot_data <- plot_data
    }
    if(!is.na(argv$prefix)){
        plot_data <- plot_data %>% filter(sample == name_df[i,1])
    }

    # EBV stat table --
    plot_data %>% group_by(cluster)%>%
        summarize(`cell with EBV` = sum(EBV_UMI>0),
                  CellNumber = length(EBV_UMI))%>%
        mutate(`cell with EBV percent` = paste(100 * round(`cell with EBV`/CellNumber,4),"%")) -> stat_df
    write.table(stat_df, paste0(outdir, "/", outfile_name, "_EBV_stat.xls"), sep="\t", row.names = F, quote = F) 
    
    # EBV status plot --
    if(length(unique(plot_data$EBV_status))==2){
        plot_color = c("red","lightgrey")
    }
    if(length(unique(plot_data$EBV_status))==1){
        plot_color = c("lightgrey")
    }
    if(argv$sort %in% c("T","True","TRUE")){
        p1 = ggplot(plot_data %>% arrange(desc(EBV_status)),
            aes(x = UMAP_1, y = UMAP_2, color = EBV_status)) +
            geom_point(size=0.1) + scale_color_manual(values = plot_color)+
            theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(outfile_name, " : EBV_status"))+
            guides(colour = guide_legend(override.aes = list(size=3)))
        ggsave(paste0(outdir, "/", outfile_name, "_EBV_status.png"), plot = p1, width = 7, height = 5)
    }else{
        p1 = ggplot(plot_data,
            aes(x = UMAP_1, y = UMAP_2, color = EBV_status)) +
            geom_point(size=0.1) + scale_color_manual(values = plot_color)+
            theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(outfile_name, " : EBV_status"))+
            guides(colour = guide_legend(override.aes = list(size=3)))
        ggsave(paste0(outdir, "/", outfile_name, "_EBV_status.png"), plot = p1, width = 7, height = 5)
    }
     
    # cell type plot --
    label_data <- plot_data %>% group_by(cluster) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    p2 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
      geom_point(size=0.1)+
      ggrepel::geom_text_repel(aes(label = cluster), data = label_data, color="black", size=3) + 
      theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(outfile_name, " : cell_type"))+
      guides(colour = guide_legend(override.aes = list(size=3)))
    ggsave(paste0(outdir, "/", outfile_name, "_cell_type.png"), plot = p2, width = 7, height = 5)

    # seurat cluster plot --
    if(!is.na(argv$anno_file)){
        label_data_1 <- plot_data %>% group_by(seurat_clusters) %>%
            summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
        plot_data$seurat_clusters <- as.character(plot_data$seurat_clusters)  # to character, for plot color
        p3 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
            geom_point(size=0.1)+
            ggrepel::geom_text_repel(aes(label = seurat_clusters), data = label_data_1, color="black", size=3) + 
            theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(outfile_name, " : seurat_clusters"))+
            guides(colour = guide_legend(override.aes = list(size=3)))
        ggsave(paste0(outdir, "/", outfile_name, "_seurat_clusters.png"), plot = p3, width = 7, height = 5)
    
        write.table(anno_file, paste0(outdir, "/", outfile_name, "_anno_file.tsv"), sep="\t", row.names =F, quote = F) 
    }
}

if(!is.na(argv$prefix)){
    name_df = data.frame(zl_name = prefix,   #  fix,  raw:unique(data_seurat$sample)
                         fj_name = fj_name)
}
if(is.na(argv$prefix)){
    name_df = data.frame(zl_name = "_",
                         fj_name = fj_name)
}
print(name_df)
out <- lapply(1:nrow(name_df),function_EBV_plot_update)


# end --
cat("EBV stat done.\n")
cat("--------\n")