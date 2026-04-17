suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds")
argv <- add_argument(argv,"--mode", help="Optional, if all, will add analysis: all_high vs all_low")
argv <- add_argument(argv,"--outdir", help="outdir, such as: sweet_tag")
argv <- parse_args(argv)

rds <- argv$rds
outdir <- argv$outdir

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = T)
}
diff_outdir <- paste0(argv$outdir, "/diff/")
if(!dir.exists(diff_outdir)){
    dir.create(diff_outdir, recursive = T)
}

#### analysis ----
data_seurat <- readRDS(rds)

# all_high vs all_low
if(!is.na(argv$mode)){
    # mkdir 
    if(!dir.exists(paste0(argv$outdir, "/diff_all/"))){
        dir.create(paste0(argv$outdir, "/diff_all/"), recursive = T)
    }
    Idents(data_seurat) <- data_seurat$sweet_tag_group
    print(table(Idents(data_seurat)))
    df <- FindAllMarkers(data_seurat, only.pos = T,verbose = F)
    df_split <- split(df, df$cluster)
    
    sapply(names(df_split),function(x){
        write.table(df_split[[x]][,c(7,1:6)], paste0(outdir, "/diff_all/",x,".tsv"), sep="\t", quote = F, row.names = F)
    })
    # diffgenes_list
    diffgenes_list <- data.frame(cluster_name = names(df_split),
                                 diffgenes_path = paste0(outdir, "/diff_all/", names(df_split), ".tsv"))
    write.table(diffgenes_list, paste0(outdir, "/diff_all/diffgenes_list.tsv"), sep="\t", row.names = F, col.names = F, quote = F)
}

# cluster_high vs cluster_low
data_seurat_list <- SplitObject(data_seurat, split.by = "cluster")
diffgenes_list <- list()

for(i in names(data_seurat_list)){
    data_analysis <- data_seurat_list[[i]]
    Idents(data_analysis) <- paste0(data_analysis$cluster, "_", data_analysis$sweet_tag_group)
    print(table(Idents(data_analysis)))
    # skip small sample
    if(min(table(Idents(data_analysis))) < 3 | length(unique(Idents(data_analysis)))==1 ) {
        print(paste0('Warning: ', i,' has not enough cells, skipping'))
        next
    }

    df <- FindAllMarkers(data_analysis, only.pos = T, verbose = F)
    df_split <- split(df, df$cluster)
    
    sapply(names(df_split),function(x){
        write.table(df_split[[x]][,c(7,1:6)], paste0(diff_outdir, "/", x, ".tsv"), sep="\t", quote = F, row.names = F)
    })

    diffgenes_list[[i]] <- data.frame(cluster_name = names(df_split),
                                      diffgenes_path = paste0(outdir, "/diff/", names(df_split), ".tsv"))
}
diffgenes_list <- do.call(rbind, diffgenes_list)
write.table(diffgenes_list, paste0(outdir, "/diff/diffgenes_list.tsv"), sep="\t", row.names = F, col.names = F, quote = F)

print("Done.")
