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
argv <- add_argument(argv, "--mod", help = "T or F")
argv <- add_argument(argv, "--subset", help = "subset")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)
cat(paste0("rds: ", argv$rds, "\n--------\n"))

if(is.na(argv$outdir)){
    argv.outdir = "./"
}
# default value -- 
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir)  
}
outdir = argv$outdir

# read rds -- 
data_seurat <- readRDS(argv$rds)
if(argv$mod == 'T'){
    data_seurat$cluster = data_seurat$cluster_singleR
}

if(!is.na(argv$subset)){
    subset_sample <- unlist(strsplit(argv$subset, split = ","))
    data_seurat <- subset(data_seurat, subset = sample %in% subset_sample)   
}

function_cell_components <- function(i){
    plot_data = data_seurat@reductions$umap@cell.embeddings %>% 
      as.data.frame() %>% 
      cbind(barcode = colnames(data_seurat))
    plot_data <- cbind(plot_data, data_seurat@meta.data)
    plot_data <- plot_data %>% filter(sample == name_df[i,1])
    
    # stat table --
    plot_data %>% group_by(cluster)%>%
        summarize(CellNumber = length(cluster))%>%
        mutate(cell_percent = paste(100 * round(CellNumber/nrow(plot_data),4),"%")) -> stat_df
    write.table(stat_df, paste0(outdir, "/", name_df[i,1], "_cell_percent_stat.xls"), sep="\t", row.names = F, quote = F) 
    
    # cell type plot --
    label_data <- plot_data %>% group_by(cluster) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    p2 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
      geom_point(size=0.1)+
      ggrepel::geom_text_repel(aes(label = cluster), data = label_data, color="black", size=3) + 
      theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(name_df[i,2], " : cell_type"))+
      guides(colour = guide_legend(override.aes = list(size=3)))
    ggsave(paste0(outdir, "/", name_df[i,1], "_cell_type.png"), plot = p2, width = 7, height = 5)
}

name_df = data.frame(zl_name = unique(data_seurat$sample),
                     zl_name_bak = unique(data_seurat$sample))
print(name_df)
out <- lapply(1:nrow(name_df),function_cell_components)
# end --
cat("cell components done.\n--------")