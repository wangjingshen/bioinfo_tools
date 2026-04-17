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
argv <- add_argument(argv, "--anno_mode", help = "cloud anno, Default: F")
argv <- add_argument(argv, "--mod", help = "same or no. Deafault: no")
argv <- add_argument(argv, "--anno_file", help = "anno file, for un-anno rds")
argv <- add_argument(argv, "--singler_anno_file", help = "singler anno file, for un-anno rds")
argv <- add_argument(argv, "--version", help = "version. v1 or singleR or v2 or rd_sr. v2 for lims. Default: v2")
#argv <- add_argument(argv, "--cele_version", help = "cele_version. Default: 1.5.1")
argv <- add_argument(argv, "--fj_path", help = "fj path, split by ,")
argv <- add_argument(argv, "--fj_name", help = "fj name, split by ,")
argv <- add_argument(argv, "--prefix", help = "prefix of cell name of zl, split by ,")
argv <- add_argument(argv, "--subset", help = "subset sample, split by ,")
argv <- add_argument(argv, "--sort", help = "plot cells in order of expression. Default:F")
#argv <- add_argument(argv, "--clean_name", help = "clean name, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: HBV_stat")
argv <- parse_args(argv)

# default value --
argv$version <- ifelse(is.na(argv$version), "v2", argv$version)
argv$mod <- ifelse(is.na(argv$mod), "no", argv$mod)
argv$order <- ifelse(is.na(argv$sort), "F", argv$sort)
argv$outdir <- ifelse(is.na(argv$outdir), "HBV_stat/", argv$outdir)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)
}
# version = argv$version # version is R function, not set, use argv$version
print(paste0("HBV stat version: ", argv$version))  # print version 

# fj_path 
if(argv$mod == "same"){
    fj_name <- unlist(strsplit(argv$fj_name, split = ","))
    fj_path <- paste0(argv$fj_path, "/", fj_name)
}
if(argv$mod == "no"){
    fj_path <- unlist(strsplit(argv$fj_path, split = ","))
    fj_name <- unlist(strsplit(argv$fj_name, split = ","))
}
print(fj_name)
print(fj_path)

if(!is.na(argv$prefix)){
    prefix = unlist(strsplit(argv$prefix, split = ","))
}


# read rds -- 
data_seurat <- readRDS(argv$rds)
if(!is.na(argv$anno_file)){
    # anno for single un-anno rds, analysis as v1  !!!
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
if(!is.na(argv$singler_anno_file)){
    anno_file <- read.table(argv$singler_anno_file, header=T,sep="\t")
    if(identical(row.names(anno_file), colnames(data_seurat))){
        data_seurat$cell_type <- anno_file$celltype
    }else{
        print("anno file barcode differernt seurat barcode")
        quit()
    }
}


if(argv$version %in% c("v1", "singleR")){
    data_seurat$sample = data_seurat$orig.ident   # add sample col
    data_seurat$cluster = data_seurat$cell_type   # add cluster col
}

if(argv$version %in% c("rd_sr")){
    data_seurat$sample <- data_seurat$orig.ident
    data_seurat$cluster <- data_seurat$celltype
}

if(!is.na(argv$subset)){
    subset_sample <- unlist(strsplit(argv$subset, split = ","))
    data_seurat <- subset(data_seurat, subset = sample %in% subset_sample)   
}

function_read_fj_csv <- function(i){
    path = fj_path[i]
    virus_tsne <- dir(path, pattern = "*tsne*", recursive = T, full.names = T)
    virus_tsne <- virus_tsne[ grep("analysis", virus_tsne)]
    if(grepl("tsv", virus_tsne)){
        fj <- read.table(virus_tsne, sep="\t", header = T, row.names = 1)
        fj$UMI[ is.na(fj$UMI)] =0
        fj$sum_UMI <- fj$UMI
    }else if(grepl("csv", virus_tsne)){
        fj <- read.table(virus_tsne, sep=",", header = T)
    }else{
        print("virus tsne file not found!")
        quit()
    }
    if(!is.na(argv$prefix)){
        bc_prefix = prefix[i]
        if(argv$version!="rd_sr"){
            fj$barcode <- paste0(bc_prefix,"_",fj$barcode)
        }
    }
    row.names(fj) <- fj$barcode
    return(fj[,c("barcode","sum_UMI")])
}
fj <- do.call(rbind,lapply(1:length(fj_name),function_read_fj_csv))
if(length(intersect(colnames(data_seurat),fj$barcode)) != ncol(data_seurat)){
    print("zl and fj not consist.")
    print(colnames(data_seurat)[1:3])
    print(fj$barcode[1:3])
    print(paste0("intersect: ", length(intersect(colnames(data_seurat),fj$barcode))))
    print(paste0("rds ncol: ", ncol(data_seurat)))
    print(paste0("fj nrow: ", nrow(fj)))
    quit()
}else{
    fj <- fj[ colnames(data_seurat),]  # match 
    data_seurat$HBV_UMI <- fj[,2]    # support 1.13.0;  raw code: data_seurat$HBV_UMI <- fj$UMI   #data_seurat$HBV_UMI <- fj$sum_UMI
    data_seurat$HBV_status <- factor(ifelse(data_seurat$HBV_UMI ==0, "HBV-", "HBV+"), levels = c("HBV+","HBV-"))
}

# HBV stat --
function_HBV_plot_update <- function(i){
    outfile_name <- ifelse(is.na(argv$prefix), name_df[i,2], name_df[i,1])
    plot_data = data_seurat@reductions$umap@cell.embeddings %>% 
      as.data.frame() %>% 
      cbind(barcode = colnames(data_seurat))
    plot_data <- cbind(plot_data, data_seurat@meta.data)
    if(!is.na(argv$prefix)){
        plot_data <- filter(plot_data, sample == name_df[i,1])
    }
    # HBV stat table --
    plot_data %>% 
        group_by(cluster) %>%
        summarize(`cell with HBV` = sum(HBV_UMI>0),
                  CellNumber = length(HBV_UMI)) %>%
        mutate(`cell with HBV percent` = paste(100 * round(`cell with HBV`/CellNumber,4),"%")) -> stat_df
    write.table(stat_df, paste0(argv$outdir, "/", outfile_name, "_HBV_stat.xls"), sep="\t", row.names = F, quote = F) 
    
    # HBV status plot --
    if(length(unique(plot_data$HBV_status))==2){
        plot_color <- c("red","lightgrey")
    }else{
        plot_color <- c("lightgrey")
    }

    if(argv$sort %in% c("T","True","TRUE")){
        p1 = ggplot(plot_data %>% arrange(desc(HBV_status)), aes(x = UMAP_1, y = UMAP_2, color = HBV_status)) +
            geom_point(size=0.1) + 
            scale_color_manual(values = plot_color) +
            ggtitle(paste0(outfile_name, " : HBV_status")) +
            theme_classic() + 
            theme(legend.title = element_blank()) + 
            guides(colour = guide_legend(override.aes = list(size=3)))
        ggsave(paste0(argv$outdir, "/", outfile_name, "_HBV_status.png"), plot = p1, width = 7, height = 5)
    }else{
        p1 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = HBV_status)) +
            geom_point(size=0.1) + 
            scale_color_manual(values = plot_color) +
            ggtitle(paste0(outfile_name, " : HBV_status")) +
            theme_classic() + theme(legend.title = element_blank()) + 
            guides(colour = guide_legend(override.aes = list(size=3)))
        ggsave(paste0(argv$outdir, "/", outfile_name, "_HBV_status.png"), plot = p1, width = 7, height = 5)
    }

    # HBV UMI plot --
    p1 = ggplot(plot_data %>% arrange(HBV_UMI), aes(x = UMAP_1, y = UMAP_2, color = HBV_UMI)) +
        geom_point(size=0.1) + 
        scale_color_gradient(low = "lightgrey", high = "red") +
        ggtitle(paste0(outfile_name, " : HBV_UMI")) +
        theme_classic() + 
        theme(legend.title = element_blank()) #+ 
        #guides(colour = guide_legend(override.aes = list(size=3)))
    ggsave(paste0(argv$outdir, "/", outfile_name, "_HBV_UMI.png"), plot = p1, width = 7, height = 5)
     
    # cell type plot --
    plot_data %>% 
        group_by(cluster) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> label_data 
    p2 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
        geom_point(size=0.1) +
        ggrepel::geom_text_repel(aes(label = cluster), data = label_data, color="black", size=3) +
        ggtitle(paste0(outfile_name, " : cell_type")) +
        theme_classic() + 
        theme(legend.title = element_blank()) + 
        guides(colour = guide_legend(override.aes = list(size=3)))
    ggsave(paste0(argv$outdir, "/", outfile_name, "_cell_type.png"), plot = p2, width = 7, height = 5)

    # seurat cluster plot --
    if(!is.na(argv$anno_file)){
        plot_data %>% group_by(seurat_clusters) %>%
            summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) -> label_data_1
        plot_data$seurat_clusters <- as.character(plot_data$seurat_clusters)  # to character, for plot color
        p3 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
            geom_point(size=0.1)+
            ggrepel::geom_text_repel(aes(label = seurat_clusters), data = label_data_1, color="black", size=3) +
            ggtitle(paste0(outfile_name, " : seurat_clusters")) +
            theme_classic() + 
            theme(legend.title = element_blank()) + 
            guides(colour = guide_legend(override.aes = list(size=3)))
        ggsave(paste0(argv$outdir, "/", outfile_name, "_seurat_clusters.png"), plot = p3, width = 7, height = 5)
        write.table(anno_file, paste0(argv$outdir, "/", outfile_name, "_anno_file.tsv"), sep="\t", row.names =F, quote = F) 
    }
}

if(is.na(argv$prefix)){
    name_df <- data.frame(zl_name = "_", fj_name = fj_name)
}else{
    name_df <- data.frame(zl_name = prefix, fj_name = fj_name)
}

#print("1")
#name_df <- ifelse(is.na(argv$prefix), data.frame(zl_name = "_", fj_name = fj_name), data.frame(zl_name = prefix, fj_name = fj_name))
#print("2")
print(name_df)
out <- lapply(1:nrow(name_df),function_HBV_plot_update)


# end --
cat("HBV stat done.\n")
cat("--------\n")