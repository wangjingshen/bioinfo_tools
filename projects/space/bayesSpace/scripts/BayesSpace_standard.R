suppressMessages({
library(BayesSpace)
library(SingleCellExperiment)
library(argparser)
library(dplyr)
library(ggplot2)
library(Seurat)
library(purrr)
library(patchwork)
library(cowplot)
library(assertthat)
library(corrplot)
library(viridis)
library(Rtsne)
})

argv <- arg_parser("")
argv <- add_argument(argv,"--datapath",help="the root of spatial data",default=NULL)
argv <- add_argument(argv,"--pca",help="Number of principal components to compute. We suggest using the top 15 PCs in most cases",default="10")
argv <- add_argument(argv,"--hvg",help="Number of highly variable genes to run PCA upon",default="2000")
argv <- add_argument(argv,"--qs",help="The values of q to evaluate",default="2,10")
argv <- add_argument(argv, "--name", help = "the sample name", default = "testsample")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)

clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
clustcol <- c(color_protocol, clustcol)
col1 <- colorRampPalette(c("#7F0000","red","red","#FF7F00","#FF7F00","yellow","yellow","cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(clustcol)
corrcol <- colorRampPalette(c("red","orange","blue","white","white"))


datapath <- argv$datapath
qs <- as.character(unlist(strsplit(argv$qs,split = ",")))
name <- argv$name

outdir <- argv$outdir
dir.create(outdir, showWarnings = FALSE)

####loading data####
sce_data <- readVisium(datapath)

# --- 核心修复 1：处理零 Count 的 spot 避免 BayesSpace 报错 ---
counts_per_spot <- Matrix::colSums(counts(sce_data))
zero_spots <- counts_per_spot == 0
if (sum(zero_spots) > 0) {
    message(paste("检测到", sum(zero_spots), "个全零 spot，正在赋予极小值以避免归一化报错..."))
    counts(sce_data)[1, zero_spots] <- 1e-6
}

seu_data <- Read10X(data.dir=file.path(datapath, "filtered_feature_bc_matrix"))

###BayesSpace###
print("################### Star  BayesSpace #############################")
###pre-processing data###
set.seed(519)
sce_data <- spatialPreprocess(sce_data, platform = 'Visium', log.normalize = TRUE, skip.PCA = FALSE, 
                              assay.type = "logcounts",
                              n.PCs = as.numeric(argv$pca),
                              n.HVGs = as.numeric(argv$hvg))

###selecting the number of clusters###
sce_data <- qTune(sce_data,qs = seq(qs[1],qs[2]),platform = 'Visium',d = as.numeric(argv$pca),
	nrep = 10000,
	burn.in = 1000
)

###find the Optima cluster number
SCE_loglik <- attr(sce_data,"q.logliks") %>% as.matrix()
diffSCE <- diff(SCE_loglik,lag=1,differences=1) %>% as.data.frame()
SCE_ck <- c(abs(diffSCE$loglik)) %>% as.data.frame()
colnames(SCE_ck)[1] <- c("k")
SCE_loglik <- as.data.frame(SCE_loglik)
SCE_loglik$sort <- rownames(SCE_loglik)
SCE_ck$sort <- rownames(SCE_ck)
SCE_loglik <- merge(SCE_loglik,SCE_ck,by="sort")
q <- SCE_loglik[which.min(SCE_loglik$k),"q"]

###spatial cluster###
set.seed(818)
bs_data <- spatialCluster(sce_data,q=q,use.dimred = "PCA",platform = 'Visium',d = as.numeric(argv$pca),init.method = 'kmeans',gamma = 2,
	model = "normal",
	nrep = 10000,
	burn.in = 1000,
	mu0 = NULL,
	lambda0 = NULL,
	alpha = 1,
	beta = 0.01,
	save.chain = FALSE,
	chain.fname = NULL
)

level_normal <- as.character(unique(bs_data@colData$spatial.cluster))[order(as.character(unique(bs_data@colData$spatial.cluster)),decreasing=F)]
levels(bs_data) <- level_normal
levels(bs_data@colData$spatial.cluster) <- level_normal

###make predict cluster for stlearn-cell_cell###
mdata <- bs_data@colData
pre_cluster <- c()
for(i in 1:length(unique(mdata$spatial.cluster))){
	pre_cluster <- c(pre_cluster,paste0("c",i))
}
tmp <- as.data.frame(matrix(NA,nrow=nrow(mdata),ncol = length(unique(mdata$spatial.cluster))))
for(i in 1:length(unique(mdata$spatial.cluster))){
tmp[i] <- ifelse(mdata$spatial.cluster==i,1,0)
}
for(i in 1:length(unique(mdata$spatial.cluster))){
	colnames(tmp)[i] = pre_cluster[i]
}
tmp$cluster_prediction <- paste0("c",mdata$spatial.cluster)
row.names(tmp) <- rownames(bs_data@colData)
write.csv(tmp,file=paste0(outdir,"/",name,".normal","_celltype_prediction.csv"), quote = FALSE,row.names = T)


####Seurat###
print("################### Star  Seurat #############################")
seu_data <- CreateSeuratObject(counts = seu_data, assay = "Spatial", project = name, min.cells = 5, names.delim="-")
image <- Read10X_Image(image.dir = file.path(datapath, "spatial"), filter.matrix = TRUE)
image <- image[Cells(x = seu_data)]
DefaultAssay(object = image) <- "Spatial"
seu_data[["slice1"]] <- image
rm(image)

mito.genes <- grep(pattern = "^MT-", x=rownames(x=seu_data[["Spatial"]]@data), value = TRUE, ignore.case = TRUE)
seu_data[["percent.mito"]] <- PercentageFeatureSet(seu_data, features = mito.genes)
Ribosomal <- grep(pattern = "^(RPL|RPS)", x = rownames(x = seu_data[["Spatial"]]@data), value = TRUE, ignore.case = TRUE)
seu_data[["percent.Ribo"]] <- PercentageFeatureSet(seu_data, features = Ribosomal)

# --- 核心修复 2：在任何画图发生前，替换 Seurat 对象中的 NA/NaN 为 0 ---
seu_data@meta.data$percent.mito[is.na(seu_data@meta.data$percent.mito)] <- 0
seu_data@meta.data$percent.Ribo[is.na(seu_data@meta.data$percent.Ribo)] <- 0

###1.QC###
qc_dir <- file.path(outdir, "1.QC")
dir.create(qc_dir, showWarnings = FALSE)
p1 <- VlnPlot(object = seu_data, features = c("percent.mito"), adjust = 1, pt.size = 0,combine = TRUE, assay = "Spatial", cols = "greenyellow") + NoLegend() + labs(x = "")
p2 <- VlnPlot(object = seu_data, features = c("percent.Ribo"), adjust = 1, pt.size = 0,combine = TRUE, assay = "Spatial", cols = "darkorange") + NoLegend() + labs(x = "")
p3 <- VlnPlot(object = seu_data, features = c("nFeature_Spatial"), adjust = 1, pt.size = 0, combine = TRUE, assay = "Spatial", cols = "darkslategray") + labs(title = "nGene", x = "") + NoLegend()
p4 <- VlnPlot(object = seu_data, features = c("nCount_Spatial"), adjust = 1, pt.size = 0, combine = TRUE, assay = "Spatial", cols = "orangered") + labs(title = "nUMI", x = "") + NoLegend()

pdf(paste0(qc_dir, '/', name, '_qc_vlnplot.pdf'), width = 6, height = 6)
print(plot_grid(p1, p2, p3, p4, ncol = 4))
dev.off()
png(paste0(qc_dir, '/', name, '_qc_vlnplot.png'), width = 6, height = 6,res=300,units='in')
print(plot_grid(p1, p2, p3, p4, ncol = 4))
dev.off()

# 优化后的分布统计函数，自动过滤 NA/NaN，且运行速度更快
mito.distrib <- function(sample, x){
  sample_clean <- sample[!is.na(sample)]
  if(length(sample_clean) == 0) return(0)
  minors <- sum(sample_clean > x)
  return (minors / length(sample_clean))
}

# 计算线粒体比例并导出
percent.mito <- Matrix::colSums(seu_data[["Spatial"]]@counts[mito.genes,]) / Matrix::colSums(seu_data[["Spatial"]]@counts)
percent.mito[is.na(percent.mito)] <- 0 
u <- as.data.frame(percent.mito)
mito <- u$percent.mito

p05 <- mito.distrib(mito, 0.1) * 100
p1 <- mito.distrib(mito, 0.2) * 100
p15 <- mito.distrib(mito, 0.3) * 100
p2 <- mito.distrib(mito, 0.5) * 100
percent_mito <- c(">10%", ">20%", ">30%", ">50%")
percent_spot_num <- c(p05, p1, p15, p2)
data.mito <- data.frame(mito_percent = percent_mito, spot_percent = percent_spot_num)
write.table(data.mito, file= paste0(qc_dir, '/', name, '_percent_mito.xls'), sep = '\t', quote = F, row.names = F)

seu_data <- NormalizeData(seu_data, assay = "Spatial", verbose = FALSE)
all_genes <- rownames(seu_data@assays$Spatial@data)
seu_data <- ScaleData(seu_data, assay = "Spatial", verbose = FALSE, features = all_genes)
seu_data  <- FindVariableFeatures(seu_data, selection.method = "vst", nfeatures = 2000)

###2.Cluster
cluster_dir <- file.path(outdir, "2.Cluster")
dir.create(cluster_dir, showWarnings = FALSE)
seu_data <- RunPCA(seu_data, assay = "Spatial", verbose = TRUE, ndims.print = 1:5, nfeatures.print = 5)

pdf(paste0(cluster_dir, '/', name, '_PCElbowPlot.pdf'))
ElbowPlot(seu_data, ndims = as.numeric(argv$pca))
dev.off()
png(paste0(cluster_dir, '/', name, '_PCElbowPlot.png'))
ElbowPlot(seu_data, ndims = as.numeric(argv$pca))
dev.off()

seu_data <- FindNeighbors(seu_data, reduction = "pca", dims = 1:as.numeric(argv$pca))
seu_data <- FindClusters(seu_data, verbose = FALSE, resolution = 0.8)
seu_data <- RunUMAP(seu_data, reduction = "pca", dims = 1:as.numeric(argv$pca), verbose = FALSE)

# --- 核心修复 3：自定义 Rtsne 去除重复坐标报错 ---
pca_emb <- seu_data@reductions$pca@cell.embeddings[, 1:as.numeric(argv$pca), drop = FALSE]
pca_emb_round <- round(pca_emb, 6)
uniq_idx <- !duplicated(pca_emb_round)
pca_unique <- pca_emb[uniq_idx, , drop = FALSE]
cells_unique <- rownames(pca_unique)

set.seed(42)
tsne_res <- Rtsne(
  pca_unique,
  dims = 2,
  pca = FALSE,
  check_duplicates = FALSE
)
tsne_mat <- matrix(NA_real_, nrow = nrow(pca_emb), ncol = 2)
rownames(tsne_mat) <- rownames(pca_emb)
colnames(tsne_mat) <- paste0("tSNE_", 1:2)
tsne_mat[cells_unique, ] <- tsne_res$Y
pca_round_str <- apply(pca_emb_round, 1, paste, collapse = "_")
first_idx <- match(pca_round_str, pca_round_str[uniq_idx])
tsne_mat[is.na(tsne_mat[,1]), ] <- tsne_res$Y[first_idx[is.na(tsne_mat[,1])], ]
tsne_dr <- CreateDimReducObject(
  embeddings = tsne_mat,
  key = "TSNE_",
  assay = DefaultAssay(seu_data)
)
seu_data@reductions$tsne <- tsne_dr

# 合并 BayesSpace 结果到 Seurat
bsclust <- data.frame(rownames(bs_data@colData))
bsclust$bscluster <- bs_data@colData$spatial.cluster
colnames(bsclust) <- c('barcode','bscluster')
spclust <- data.frame(rownames(seu_data@meta.data))
colnames(spclust) <- c('barcode')
spclust <- left_join(spclust,bsclust,by='barcode')
seu_data@meta.data$bscluster <- spclust$bscluster
seu_data <- SetIdent(seu_data,value=seu_data@meta.data$bscluster)
levels(seu_data) <- sort(levels(seu_data))

# 动态保护颜色数量（确保不会因为 Cluster 太多而报错）
num_clusters <- length(levels(seu_data))
if(length(clustcol) < num_clusters) {
  clustcol <- colorRampPalette(clustcol)(num_clusters)
}

#####保存RDS文件
rds_dir <- file.path(outdir, "0.Rds")
dir.create(rds_dir, showWarnings = FALSE)
saveRDS(seu_data, file = paste0(rds_dir, "/", name, ".rds"))

clust <- summary(seu_data@active.ident)
cluster_spot <- as.data.frame(clust)
spot_number <- sum(cluster_spot$clust)
print(paste0("Total spot number used: ", spot_number))

pt_use <- 0.6
if(spot_number > 1000) pt_use <- 0.4
if(spot_number > 2500) pt_use <- 0.3
if(spot_number > 4000) pt_use <- 0.2
if(spot_number > 5500) pt_use <- 0.15
if(spot_number > 6500) pt_use <- 0.1

p1 <- FeaturePlot(seu_data, reduction = "tsne", features = c("percent.mito"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q95")
p2 <- FeaturePlot(seu_data, reduction = "tsne", features = c("percent.Ribo"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99")
p3 <- FeaturePlot(seu_data, reduction = "tsne", features = c("nFeature_Spatial"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99") + labs(title = "nGene")
p4 <- FeaturePlot(seu_data, reduction = "tsne", features = c("nCount_Spatial"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99") + labs(title = "nUMI")
pdf(paste0(qc_dir, '/', name, '_qc_tsneFeatureplot.pdf'), width = 6, height = 6)
print(plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()
png(paste0(qc_dir, '/', name, '_qc_tsneFeatureplot.png'), width = 6, height = 6,res=300,units='in')
print(plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()

p1 <- FeaturePlot(seu_data, reduction = "umap", features = c("percent.mito"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99")
p2 <- FeaturePlot(seu_data, reduction = "umap", features = c("percent.Ribo"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99")
p3 <- FeaturePlot(seu_data, reduction = "umap", features = c("nFeature_Spatial"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99") + labs(title = "nGene")
p4 <- FeaturePlot(seu_data, reduction = "umap", features = c("nCount_Spatial"), cols = rev(col1(9)), pt.size = pt_use/2, min.cutoff = "q1", max.cutoff = "q99") + labs(title = "nUMI")
pdf(paste0(qc_dir, '/', name, '_qc_umapFeatureplot.pdf'), width = 6, height = 6)
print(plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()
png(paste0(qc_dir, '/', name, '_qc_umapFeatureplot.png'), width = 6, height = 6,res=300,units='in')
print(plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()

p1 <- SpatialFeaturePlot(seu_data, features = c("percent.mito"), min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "percent.mito")
p2 <- SpatialFeaturePlot(seu_data, features = c("percent.Ribo"), min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "percent.Ribo")
p3 <- SpatialFeaturePlot(seu_data, features = c("nFeature_Spatial"), min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "nGene")
p4 <- SpatialFeaturePlot(seu_data, features = c("nCount_Spatial"), min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "nUMI")
pdf(paste0(qc_dir, '/', name, '_qc_Tissueplot.pdf'), width = 6, height = 6)
print(plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()
png(paste0(qc_dir, '/', name, '_qc_Tissueplot.png'), width = 6, height = 6,res=300,units='in')
print(plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()

ratio <- sprintf('%.2f', 100 * cluster_spot$clust / sum(cluster_spot$clust))
regions <- paste0('cluster ', levels(x = seu_data))
labels <- paste0(regions, '(', ratio, '%)')
data <- data.frame(values = ratio, type = regions)
data$type <- factor(data$type, levels = (regions))
rl <- unique(ratio)
data$values <- factor(data$values, levels = (rl))

p <- ggplot(data,aes(x = '', y=values, fill = type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = 'y', start = 0) +
  labs(x = '', y = '', fill = 'type', title = paste0('Percent of spots number per Cluster ', '(', name, ')')) +
  scale_fill_manual(name = '', values = clustcol, limits = rev(regions), labels = labels) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom', legend.direction = 'vertical',
        panel.background = element_rect(fill = 'transparent', colour = NA),
        axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pdf(paste0(cluster_dir, '/', name, '_PercentPerSpot.pdf'))
print(p)
dev.off()
png(paste0(cluster_dir, '/', name, '_PercentPerSpot.png'))
print(p)
dev.off()

cluster.averages0 <- AverageExpression(object = seu_data, return.seurat = TRUE, verbose = FALSE)
cluster.averages1 <- AverageExpression(object = seu_data)
cluster.averages <- cluster.averages1$Spatial
colnames(cluster.averages) <- regions
write.table(cluster.averages, file= paste0(cluster_dir, '/', name, '_cluster_averages.xls'), sep = '\t', quote = FALSE, row.names = TRUE)

#####相关性
M <- cor(cluster.averages)
order.hc <- corrMatOrder(M, order = "hclust")
M.hc <- M[order.hc, order.hc]
pdf(paste0(cluster_dir, '/', name, '_corrplot.pdf'))
corrplot(M.hc, method = "number", cl.lim = c(0, 1), tl.col = "black", col = rev(corrcol(50)))
dev.off()
png(paste0(cluster_dir, '/', name, '_corrplot.png'))
corrplot(M.hc, method = "number", cl.lim = c(0, 1), tl.col = "black", col = rev(corrcol(50)))
dev.off()

#####降维图
spot_type <- as.data.frame(seu_data@active.ident)
colnames(spot_type) <- c("spot_type")
write.table(spot_type, file = paste0(cluster_dir, '/', name, '_spot_type.xls'), sep = '\t', quote = FALSE, row.names = TRUE)

p <- DimPlot(seu_data, reduction = "pca",group.by='bscluster', pt.size = pt_use, cols = clustcol) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_pca.pdf'))
print(p)
dev.off()
png(paste0(cluster_dir, '/', name, '_pca.png'))
print(p)
dev.off()

p <- DimPlot(seu_data, reduction = "tsne",group.by='bscluster', pt.size = pt_use, cols = clustcol) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_tsne.pdf'))
print(p)
dev.off()
png(paste0(cluster_dir, '/', name, '_tsne.png'))
print(p)
dev.off()

p <- DimPlot(seu_data, reduction = "umap",group.by='bscluster', pt.size = pt_use, cols = clustcol) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_umap.pdf'))
print(p)
dev.off()
png(paste0(cluster_dir, '/', name, '_umap.png'))
print(p)
dev.off()

cluster_order = levels(seu_data@meta.data$bscluster)
color_order = clustcol[1:length(cluster_order)]
cols <- list()
for(i in 1:length(cluster_order)) {
  cols[[cluster_order[i]]] <- color_order[i]
}
p <- SpatialDimPlot(seu_data,group.by='bscluster', cols = cols) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_tissue.pdf'))
print(p)
dev.off()
png(paste0(cluster_dir, '/', name, '_tissue.png'))
print(p)
dev.off()

p_elbow <- qPlot(sce_data)
pdf(paste0(cluster_dir,'/',name,'.elbow_cluster.pdf'),width = 7,height = 7)
print(p_elbow)
dev.off()
png(paste0(cluster_dir,'/',name,'.elbow_cluster.png'),width = 7,height = 7,res = 300,units = 'in')
print(p_elbow)
dev.off()


#####寻找差异基因####
diff_dir <- file.path(outdir, "3.Diff")
dir.create(diff_dir, showWarnings = FALSE)
diffgene_dir <- file.path(diff_dir, "3.1.diffgene")
dir.create(diffgene_dir, showWarnings = FALSE)

# --- 核心修复 4：正确遍历实际存在的 cluster 名字 ---
all_clusters <- levels(seu_data)
for (clust_id in all_clusters) {
  cluster.markers <- FindMarkers(
    object = seu_data,
    ident.1 = clust_id,
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  if ("avg_log2FC" %in% colnames(cluster.markers)) {
    colnames(cluster.markers)[colnames(cluster.markers) == "avg_log2FC"] <- "avg_logFC"
  }
  
  cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE), ]
  write.table(
    cluster.markers,
    file = file.path(diffgene_dir, paste0(name, '_cluster_', clust_id, '_diffgenes.xls')),
    sep = '\t',
    quote = FALSE,
    row.names = TRUE
  )
}

seu_data.marker <- FindAllMarkers(object = seu_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
if ("avg_log2FC" %in% colnames(seu_data.marker)) {
    colnames(seu_data.marker)[colnames(seu_data.marker) == "avg_log2FC"] <- "avg_logFC"
}
markergene <- seu_data.marker %>% group_by(cluster) %>% top_n(2, avg_logFC)
markergenetop10 <- seu_data.marker %>% group_by(cluster) %>% top_n(10, avg_logFC)

diffheatmap_dir <- file.path(diff_dir, "3.2.heatmap")
dir.create(diffheatmap_dir, showWarnings = FALSE)

pdf(paste0(diffheatmap_dir, '/', name, '_diffgene_DOHeatmapplot.pdf'))
print(DoHeatmap(seu_data, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis())
dev.off()
png(paste0(diffheatmap_dir, '/', name, '_diffgene_DOHeatmapplot.png'))
print(DoHeatmap(seu_data, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis())
dev.off()

pdf(paste0(diffheatmap_dir, '/', name, '_diffgenetop10_DOHeatmapplot.pdf'))
print(DoHeatmap(seu_data, features = markergenetop10$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 2)) + scale_fill_viridis())
dev.off()
png(paste0(diffheatmap_dir, '/', name, '_diffgenetop10_DOHeatmapplot.png'))
print(DoHeatmap(seu_data, features = markergenetop10$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 2)) + scale_fill_viridis())
dev.off()

Allmarker <- seu_data.marker[order(seu_data.marker$avg_logFC, decreasing = TRUE),]
Allmarker <- Allmarker[,!colnames(Allmarker) %in% c("gene")]
write.table(Allmarker, file = paste0(diffgene_dir, '/', name, '_all_diffgenes.xls'), sep = '\t', quote = FALSE, row.names = TRUE)

allmarker <- as.data.frame(markergene)
marker <- allmarker$gene
marker <- unique(marker)
print(paste0("Total Unique Markers used for Dot/Vln plots: ", length(marker)))

pdf(paste0(diffheatmap_dir, '/', name, '_cluster_Heatmapplot.pdf'))
print(DoHeatmap(cluster.averages0, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = FALSE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis())
dev.off()
png(paste0(diffheatmap_dir, '/', name, '_cluster_Heatmapplot.png'))
print(DoHeatmap(cluster.averages0, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = FALSE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis())
dev.off()

######差异基因DotPlot
diffdotplot_dir <- file.path(diff_dir, "3.3.dotplot")
dir.create(diffdotplot_dir, showWarnings = FALSE)
x_size <- 9
dot_size <- 7
if (length(marker) > 10) dot_size <- 6
if (length(marker) > 15) dot_size <- 6
if (length(marker) > 20) dot_size <- 5
if (length(marker) > 25) dot_size <- 4
if (length(marker) > 30) dot_size <- 3
x_size <- dot_size + 2.2
lab <- length(levels(seu_data))
l_size <- 10
if (lab > 10) l_size <- 9
if (lab > 15) l_size <- 8
if (lab > 20) l_size <- 7
if (lab > 25) l_size <- 6
if (lab > 30) l_size <- 5

pdf(paste0(diffdotplot_dir, '/', name, '_diffgene_Dotplot.pdf'))
print(DotPlot(object = seu_data, features = marker, cols = c("blue", "red"), assay = "Spatial", dot.scale = dot_size) + RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(color = guide_colorbar(title = "avg.exp.scale"), size = guide_legend(title = "pct.exp")) + theme(axis.text.x = element_text(size = x_size), axis.text.y = element_text(size = l_size)))
dev.off()
png(paste0(diffdotplot_dir, '/', name, '_diffgene_Dotplot.png'))
print(DotPlot(object = seu_data, features = marker, cols = c("blue", "red"), assay = "Spatial", dot.scale = dot_size) + RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(color = guide_colorbar(title = "avg.exp.scale"), size = guide_legend(title = "pct.exp")) + theme(axis.text.x = element_text(size = x_size), axis.text.y = element_text(size = l_size)))
dev.off()

###差异基因VlnPlot
diffvlnplot_dir <- file.path(diff_dir, "3.4.violinplot")
dir.create(diffvlnplot_dir, showWarnings = FALSE)
g <- 0
pnum <- 1
plot_g <- c()
for (i  in marker){
  if(g == 4){
    p1 <- VlnPlot(object = seu_data, features = plot_g[1], pt.size = 0, cols = clustcol, ncol = 2, combine = TRUE) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()
    p2 <- VlnPlot(object = seu_data, features = plot_g[2], pt.size = 0, cols = clustcol, ncol = 2, combine = TRUE) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()
    p3 <- VlnPlot(object = seu_data, features = plot_g[3], pt.size = 0, cols = clustcol, ncol = 2, combine = TRUE) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()
    p4 <- VlnPlot(object = seu_data, features = plot_g[4], pt.size = 0, cols = clustcol, ncol = 2, combine = TRUE) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()
    pdf(paste0(diffvlnplot_dir, '/', name, '_diffgene_vlnplot', pnum, '.pdf'))
    print(plot_grid(p1, p2, p3, p4, ncol = 2))
    dev.off()
    png(paste0(diffvlnplot_dir, '/', name, '_diffgene_vlnplot', pnum, '.png'))
    print(plot_grid(p1, p2, p3, p4, ncol = 2))
    dev.off()
    pnum <- pnum + 1
    g <- 1
    plot_g <- c(i)
  }else{
    plot_g <- c(plot_g, i)
    g <- g + 1
  }
}
if (g > 0){
  plot <- list()
  for (i in 1:g){
    plot[[i]] <- VlnPlot(object = seu_data, features = plot_g[i], pt.size = 0, cols = clustcol, ncol = 2, combine = TRUE) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()
  }
  pdf(paste0(diffvlnplot_dir, '/', name, '_diffgene_vlnplot', pnum, '.pdf'))
  print(plot_grid(plotlist = plot, ncol = 2))
  dev.off()
  png(paste0(diffvlnplot_dir, '/', name, '_diffgene_vlnplot', pnum, '.png'))
  print(plot_grid(plotlist = plot, ncol = 2))
  dev.off()
  g <- 0
  plot_g <- c()
}

####差异基因FeaturePlot
difffeatureplot_dir <- file.path(diff_dir, "3.5.featureplot")
dir.create(difffeatureplot_dir, showWarnings = FALSE)
g <- 0
plot_g <- c()
pnum <- 1
for (i in marker){
  if (g==4){
    plot_tsne <- list()
    plot_umap <- list()
    plot_tissue <- list()
    for (j in 1:g){
      plot_tsne[[j]] <- FeaturePlot(object = seu_data, features = plot_g[j], cols = c("lightgrey", "red"), pt.size = pt_use/2, reduction = "tsne", min.cutoff = "q1", max.cutoff = "q99") + NoLegend()
      plot_umap[[j]] <- FeaturePlot(object = seu_data, features = plot_g[j], cols = c("lightgrey", "red"), pt.size = pt_use/2, reduction = "umap", min.cutoff = "q1", max.cutoff = "q99") + NoLegend()
      plot_tissue[[j]] <- SpatialFeaturePlot(seu_data, features = plot_g[j], min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = plot_g[j])
    }
    pdf(paste0(difffeatureplot_dir, '/', name, '_diffgene_tsneFeatureplot', pnum, '.pdf'))
    print(plot_grid(plotlist = plot_tsne, ncol = 2))
    dev.off()
    png(paste0(difffeatureplot_dir, '/', name, '_diffgene_tsneFeatureplot', pnum, '.png'))
    print(plot_grid(plotlist = plot_tsne, ncol = 2))
    dev.off()
    pdf(paste0(difffeatureplot_dir, '/', name, '_diffgene_umapFeatureplot', pnum, '.pdf'))
    print(plot_grid(plotlist = plot_umap, ncol = 2))
    dev.off()
    png(paste0(difffeatureplot_dir, '/', name, '_diffgene_umapFeatureplot', pnum, '.png'))
    print(plot_grid(plotlist = plot_umap, ncol = 2))
    dev.off()
    pdf(paste0(difffeatureplot_dir, '/', name, '_diffgene_Tissueplot', pnum, '.pdf'))
    print(plot_grid(plotlist = plot_tissue, ncol = 2))
    dev.off()
    png(paste0(difffeatureplot_dir, '/', name, '_diffgene_Tissueplot', pnum, '.png'))
    print(plot_grid(plotlist = plot_tissue, ncol = 2))
    dev.off()
    pnum <- pnum + 1
    g <- 1
    plot_g <- c(i)
  }else{
    plot_g <- c(plot_g,i)
    g <- g + 1
  }
}
if(g > 0){
  plot_tsne <- list()
  plot_umap <- list()
  plot_tissue <- list()
  for (j in 1:g){
    plot_tsne[[j]] <- FeaturePlot(object = seu_data, features = plot_g[j], cols = c("lightgrey", "red"), pt.size = pt_use/2, reduction = "tsne", min.cutoff = "q1", max.cutoff = "q99") + NoLegend()
    plot_umap[[j]] <- FeaturePlot(object = seu_data, features = plot_g[j], cols = c("lightgrey", "red"), pt.size = pt_use/2, reduction = "umap", min.cutoff = "q1", max.cutoff = "q99") + NoLegend()
    plot_tissue[[j]] <- SpatialFeaturePlot(seu_data, features = plot_g[j], min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = plot_g[j])
  }
  pdf(paste0(difffeatureplot_dir, '/', name, '_diffgene_tsneFeatureplot', pnum, '.pdf'))
  print(plot_grid(plotlist = plot_tsne, ncol = 2))
  dev.off()
  png(paste0(difffeatureplot_dir, '/', name, '_diffgene_tsneFeatureplot', pnum, '.png'))
  print(plot_grid(plotlist = plot_tsne, ncol = 2))
  dev.off()
  pdf(paste0(difffeatureplot_dir, '/', name, '_diffgene_umapFeatureplot', pnum, '.pdf'))
  print(plot_grid(plotlist = plot_umap, ncol = 2))
  dev.off()
  png(paste0(difffeatureplot_dir, '/', name, '_diffgene_umapFeatureplot', pnum, '.png'))
  print(plot_grid(plotlist = plot_umap, ncol = 2))
  dev.off()
  pdf(paste0(difffeatureplot_dir, '/', name, '_diffgene_Tissueplot', pnum, '.pdf'))
  print(plot_grid(plotlist = plot_tissue, ncol = 2))
  dev.off()
  png(paste0(difffeatureplot_dir, '/', name, '_diffgene_Tissueplot', pnum, '.png'))
  print(plot_grid(plotlist = plot_tissue, ncol = 2))
  dev.off()
}

####空间表达
space_dir <- file.path(outdir, "5.Space")
dir.create(space_dir, showWarnings = FALSE)
seu_data <- FindSpatiallyVariableFeatures(seu_data, assay = "Spatial", features = VariableFeatures(seu_data)[1:1000], selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(seu_data, selection.method = "markvariogram"), 12)
write.table(top.features, file = paste0(space_dir, "/", name, "_SpatiallyVariableGene.txt"), sep = "\t", quote =  FALSE, row.names = FALSE, col.names = FALSE)

g <- 0
plot_g <- c()
pnum <- 1
for (i in top.features){
  if (g==4){
    plot_tissue <- list()
    for (j in 1:g){
        plot_tissue[[j]] <- SpatialFeaturePlot(seu_data, features = plot_g[j], min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = plot_g[j])
    }
    pdf(paste0(space_dir, '/', name, '_spacegene_Tissueplot', pnum, '.pdf'))
    print(plot_grid(plotlist = plot_tissue, ncol = 2))
    dev.off()
    png(paste0(space_dir, '/', name, '_spacegene_Tissueplot', pnum, '.png'))
    print(plot_grid(plotlist = plot_tissue, ncol = 2))
    dev.off()
    pnum <- pnum + 1
    g <- 1
    plot_g <- c(i)
  }else{
    plot_g <- c(plot_g, i)
    g <- g + 1
  }
}
if(g > 0){
  plot_tissue <- list()
  for (j in 1:g){
    plot_tissue[[j]] <- SpatialFeaturePlot(seu_data, features = plot_g[j], min.cutoff = "q1", max.cutoff = "q99", alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = plot_g[j])
  }
  pdf(paste0(space_dir, '/', name, '_spacegene_Tissueplot', pnum, '.pdf'))
  print(plot_grid(plotlist = plot_tissue, ncol = 2))
  dev.off()
  png(paste0(space_dir, '/', name, '_spacegene_Tissueplot', pnum, '.png'))
  print(plot_grid(plotlist = plot_tissue, ncol = 2))
  dev.off()
}

print("Workflow Completed Successfully!")