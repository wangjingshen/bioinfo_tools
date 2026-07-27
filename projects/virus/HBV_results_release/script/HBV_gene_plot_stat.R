options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--prefix", help = "prefix, for annotated rds")
argv <- add_argument(argv, "--virus_tsne", help = "virus tsne")
argv <- add_argument(argv, "--virus_matrix", help = "virus matrix")
argv <- add_argument(argv, "--HBV_positive_stat", help = "HBV positive stat")
argv <- add_argument(argv, "--reduction", help = "reduction, Default: umap")
#argv <- add_argument(argv, "--step", help = "step 1 for HBV status plot, step2 for HBV virus genes plot, step3 for stat. Default: 123")
argv <- add_argument(argv, "--outdir_02", help = "outdir 02.")
argv <- add_argument(argv, "--outdir_03", help = "outdir 03.")
argv <- add_argument(argv, "--outdir_04", help = "outdir 04.")
argv <- add_argument(argv, "--name", help = "name.")

argv <- parse_args(argv)
argv$reduction <- ifelse(is.na(argv$reduction), 'umap', argv$reduction)
#argv$step <- ifelse(is.na(argv$step), '12', argv$step)

#if(!dir.exists(argv$outdir_02)){
#    dir.create(argv$outdir_02)
#}

# start --
cat(paste0(argv$name, " starting...", "\n"))
cat("----------------\n")

if(!is.na(argv$rds)){   # argv$rds != "None", for python run R
    print(paste0("rds: ", argv$rds))
    data_seurat <- readRDS(argv$rds)
    if(!is.na(argv$prefix)){
        data_seurat <- subset(data_seurat, subset= sample == argv$prefix)  # subset
    }
}
print(paste0("virus_tsne : ", argv$virus_tsne))  # cat(paste0("virus_tsne : ", argv$virus_tsne, "\n"))
print(paste0("virus matrix : ", argv$virus_matrix))
print(paste0("celescope HBV positive stat : ", argv$HBV_positive_stat))
cat("----------------\n")

# read virus tsne --
virus_tsne <- read.table(argv$virus_tsne, sep="\t", header = T, row.names = 1)
virus_tsne$UMI[ is.na(virus_tsne$UMI)] =0
if(!is.na(argv$prefix)){
    virus_tsne$barcode <- paste0(argv$prefix, "_", virus_tsne$barcode)
}
row.names(virus_tsne) <- virus_tsne$barcode
virus_tsne <- virus_tsne[, c("barcode", "UMI")]

if(length(intersect(colnames(data_seurat), virus_tsne$barcode)) == 0){
    print("zl and fj not consist.")
    quit()
}else{
    virus_tsne <- virus_tsne[ colnames(data_seurat),]  # match 
    data_seurat$HBV_UMI <- virus_tsne[,2]
    data_seurat$HBV_status <- ifelse(data_seurat$HBV_UMI ==0, "HBV-", "HBV+")
    data_seurat$HBV_status <- factor(data_seurat$HBV_status,levels = c("HBV+", "HBV-"))
}

# virus matrix --
virus_matrix <- Read10X(argv$virus_matrix)
virus_matrix = merge(virus_matrix, data.frame(row.names = c("forS","forX","forcccDNA","forpg","forrcDNA")), by=0, all.y = T)
row.names(virus_matrix) <- virus_matrix[,1]
virus_matrix <- virus_matrix[,-1]
if(!is.na(argv$prefix)){
    colnames(virus_matrix) <- paste0(argv$prefix, "_", colnames(virus_matrix))
}
# merge barcode
virus_matrix_raw <- t(virus_matrix)   # virus_matrix_raw for 03 gene positive stat
virus_matrix_raw[is.na(virus_matrix_raw)] = 0
virus_matrix <- merge(data.frame(row.names = colnames(data_seurat)), t(virus_matrix), by = 0, all.x = T)

row.names(virus_matrix) <- virus_matrix[,1]
virus_matrix <- virus_matrix[,-1]
virus_matrix[is.na(virus_matrix)] = 0

data_seurat$forS <- virus_matrix$forS
data_seurat$forX <- virus_matrix$forX
data_seurat$forcccDNA <- virus_matrix$forcccDNA
data_seurat$forpg <- virus_matrix$forpg
data_seurat$forrcDNA <- virus_matrix$forrcDNA

#
virus_stat <- read.table(argv$HBV_positive_stat, sep=":")

# get plot data --
if(argv$reduction=="umap"){
    plot_data = data_seurat@reductions$umap@cell.embeddings %>% 
        as.data.frame() %>% 
        cbind(barcode = colnames(data_seurat))
}
if(argv$reduction=="tsne"){
    plot_data = data_seurat@reductions$tsne@cell.embeddings %>% 
        as.data.frame() %>% 
        cbind(barcode = colnames(data_seurat))
}
plot_data <- cbind(plot_data, data_seurat@meta.data)

# 02/xxx_HBV_status.pdf ----
#if(grepl('1', argv$step)){
if(length(unique(data_seurat$HBV_status))==2){
    plot_color = c("red","lightgrey")
}
if(length(unique(data_seurat$HBV_status))==1){
    plot_color = c("lightgrey")
}
p1 <- ggplot(plot_data %>% arrange(desc(HBV_status)), aes(x = UMAP_1, y = UMAP_2, color = HBV_status)) +
        geom_point(size=0.1) + scale_color_manual(values = plot_color)+
        theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(argv$name, " : HBV_status"))+
        guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(paste0(argv$outdir_02, "/", argv$name, "_HBV_status.pdf"), plot = p1, height = 5, width = 7)
#}

# 02/xxx_HBV_genes_stat.tsv ----
HBV_genes_stat <- data_seurat@meta.data %>%
    group_by(cluster) %>%
    summarise(S_sum = sum(forS),
              X_sum = sum(forX),
              cccDNA_sum = sum(forcccDNA),
              pgRNA_sum = sum(forpg),
              rcDNA_sum = sum(forrcDNA))
write.table(HBV_genes_stat, paste0(argv$outdir_02, "/", argv$name, "_HBV_genes_stat.tsv"), sep="\t", quote = F, row.names = F)

# 03/xxx_HBV_virus_genes.pdf ----
p2 <- FeaturePlot(data_seurat,cols = c("lightgrey","red"), features = c("forS","forX","forcccDNA","forpg","forrcDNA"), order = T, reduction = argv$reduction, keep.scale = 'all', pt.size = 0.1, ncol=3)
ggsave(paste0(argv$outdir_03, "/", argv$name, "_HBV_virus_genes.pdf"), plot = p2, height = 7, width = 12)

# 03/xxx_HBV_and_virus_genes_positive_stat.tsv ----
HBV_gene_positive_stat <- data.frame(
    Name = paste0("Number of cells with ", c('virus', gsub("pg", "pgRNA", gsub("for", "", names(colSums(virus_matrix_raw > 0)))))),
    CellNumber = c(as.numeric(word(gsub(",","",virus_stat[3,2]), 1, sep = fixed("("))), as.numeric(colSums(virus_matrix_raw > 0))))
HBV_gene_positive_stat$Percent <- paste0(round(HBV_gene_positive_stat$CellNumber/as.numeric(gsub(",", "", virus_stat[1,2]))*100, 3), "%")
#Percent = paste0(round(c(as.numeric(word(virus_stat[3,2], 1, sep=fixed("("))), c(colSums(virus_matrix>0)))/as.numeric(gsub(",", "", virus_stat[1,2]))*100, 3), "%"))
write.table(HBV_gene_positive_stat, paste0(argv$outdir_03, "/", argv$name, "_HBV_and_virus_genes_positive_stat.tsv"), sep="\t", quote = F, row.names = F)

# 04/xxx_HBV_stat.tsv ----
data_seurat@meta.data %>% 
    group_by(cluster)%>%
    summarize(`cell with HBV` = sum(HBV_UMI > 0),
              CellNumber = length(HBV_UMI)) %>%
    mutate(`cell with HBV percent` = paste(100 * round(`cell with HBV`/CellNumber,4),"%")) -> stat_df
write.table(stat_df, paste0(argv$outdir_04, "/", argv$name, "_HBV_cell_type_positive_stat.tsv"), sep="\t", row.names = F, quote = F) 

# 04/xxx_cell_type.pdf ----
#p3 <- DimPlot(data_seurat, group.by = "cluster", reduction = argv$reduction, pt.size = 0.1)
#ggsave(paste0(argv$outdir_04, "/", argv$name, "_cell_type.pdf"), plot = p3, height = 5, width = 7)

label_data <- plot_data %>% group_by(cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
p3 = ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size=0.1)+
    ggrepel::geom_text_repel(aes(label = cluster), data = label_data, color="black", size=3) + 
    theme_classic() + theme(legend.title = element_blank()) + ggtitle(paste0(argv$name, " : cell_type"))+
    guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(paste0(argv$outdir_04, "/", argv$name, "_cell_type.pdf"), plot = p3, width = 7, height = 5)

cat("----------------\n")
cat(paste0(argv$name, " done.", "\n"))
cat("----------------\n\n")