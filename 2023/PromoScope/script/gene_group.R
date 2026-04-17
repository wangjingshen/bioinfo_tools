suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101","OrangeRed","SlateBlue","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink","Red")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds")
argv <- add_argument(argv,"--cluster", help="cluster")
argv <- add_argument(argv,"--gene", help="gene")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- parse_args(argv)

rds <- argv$rds
cluster <- unlist(str_split(argv$cluster, ','))
gene <- unlist(str_split(argv$gene, ','))
outdir <- argv$outdir

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

data_seurat <- readRDS(argv$rds)

gene <- intersect(gene, row.names(data_seurat@assays$RNA@counts))


function_analysis <- function(data, sub_cluster, gene){
    data_sub = subset(data, subset = cluster == sub_cluster)
    data_sub$group <- as.character(data_sub$group)
    data_sub$group <- ifelse(data_sub@assays$RNA@data[ gene,] > 0, str_glue(gene,"_positive"), str_glue(gene,"_negative"))
    Idents(data_sub) <- data_sub$group
    print(sub_cluster)
    print(table(data_sub$group))
    if(length(unique(data_sub$group))==2 & min(table(data_sub$group)) >= 3){
        options(repr.plot.height= 5, repr.plot.width = 6)
        data_sub@meta.data %>%
            ggplot(aes(x = cluster, y = sweet_tag_CLR, fill = group))+
            geom_boxplot(outlier.shape = NA)+
            guides(fill = guide_legend(title = "group"))+
            scale_fill_manual(values = color_protocol)+
            ylab("sweet_tag")+
            ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test")+
            theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1 ),
                                    plot.title = element_blank(),
                                    axis.title.x = element_blank())
        ggsave(str_glue(outdir, "/BoxPlot_", sub_cluster,"_", gene,".pdf"), height = 5, width = 6)
        ggsave(str_glue(outdir, "/BoxPlot_", sub_cluster,"_", gene,".png"), height = 5, width = 6)

        saveRDS(data_sub, str_glue(outdir, "/data_", sub_cluster,"_", gene,".rds"))
    }else{
        print("skip")
    }
}

for (i in cluster){
    for (j in gene){
        function_analysis(data_seurat, i, j)
    }
}

print("Done.")