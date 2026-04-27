suppressWarnings(suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
    library(patchwork)}))

color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8",
"#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F",
"#000101","OrangeRed","SlateBlue","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink","Red","#4682B4",
"#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00",
"#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE",
"#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD",
"#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF",
"#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4",
"#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F",
"#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347",
"#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF",
"#1E90FF","#191970","#E8E8E8","#FFDAB9")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--anno_file", help = "anno file, for un-anno rds")
argv <- add_argument(argv, "--outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

# read rds -- 
data_seurat <- readRDS(argv$rds)
data_seurat$seurat_clusters <- as.character(data_seurat$seurat_clusters)

anno_file <- read.table(argv$anno_file, header=T,sep="\t")
#anno_file$cluster <- anno_file$cluster - 1    for old version
anno_file <- anno_file[!duplicated(anno_file$cluster),]    # if duplicated, use first anno
anno_file$cluster <- as.character(anno_file$cluster)
anno_file$cell_type[is.na(anno_file$cell_type)] = "NA"     # to character, for plot label "NA"; if NA, not labeled
#data_seurat@meta.data%>%
#    left_join(anno_file,by = c("seurat_clusters" = "cluster")) -> anno_df

data_seurat$cluster_auto_assign <- plyr::mapvalues(x = data_seurat$seurat_clusters, 
    from = anno_file$cluster , to = anno_file$cell_type)
p1 = DimPlot(data_seurat, group.by = "cluster_auto_assign", cols = color_protocol, label = T, repel=T)
ggsave(str_glue("{outdir}/cluster_auto_assign.png"), plot=p1, height=9, width=12)
ggsave(str_glue("{outdir}/cluster_auto_assign.pdf"), plot=p1, height=9, width=12)

# filter species if need
data_seurat$cluster <- gsub("_mmu", "", gsub("_hs", "", data_seurat$cluster_auto_assign))
p2 = DimPlot(data_seurat, group.by = "cluster", cols = color_protocol)
ggsave(str_glue("{outdir}/cluster.png"), plot=p2, height=9, width=12)
ggsave(str_glue("{outdir}/cluster.pdf"), plot=p2, height=9, width=12)

#data_seurat@meta.data %>% group_by(cluster_auto_assign)%>%
#    summarize(CellNumber = length(cluster_auto_assign))%>%
#    mutate(cell_percent = paste(100 * round(CellNumber/ncol(data_seurat),4),"%")) -> stat_df
#write.table(stat_df, paste0(outdir, "/cell_number_stat.xls"), sep="\t", row.names = F, quote = F) 


# cellular_composition
#p2 <- data_seurat@meta.data %>%
#    ggplot(aes(x=sample,fill=cluster_auto_assign))+
#        geom_bar(color="black",position="fill",width = 0.5)+
#        theme_bw()+
#        theme(legend.title = element_blank(),axis.title = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
#ggsave(paste0(outdir, "/cellular_composition.pdf"), plot= p2, width = 4, height = 5)
#ggsave(paste0(outdir, "/cellular_composition.png"), plot= p2, width = 4, height = 5)

saveRDS(data_seurat, paste0(outdir, "/seurat_anno.rds"))

# end --
cat("auto assign plot done.\n")
cat("--------\n")