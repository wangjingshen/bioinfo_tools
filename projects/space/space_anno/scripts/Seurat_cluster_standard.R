#conda activate Seuratv5
suppressMessages({
#library(BayesSpace)
#library(SingleCellExperiment)
library(argparser) #需要安装
library(dplyr)
library(ggplot2)
library(Seurat)
library(purrr)
library(patchwork)
library(cowplot)
library(assertthat)
library(stringr)
#library(corrplot)  
library(viridis) #需要安装
})

argv <- arg_parser("")
argv <- add_argument(argv,"--h5",help="spatial h5 file",default=NULL)
argv <- add_argument(argv,"--spatial",help="spatial folder",default=NULL)
argv <- add_argument(argv,"--pca",help="Number of principal components to compute. We suggest using the top 15 PCs in most cases",default="20")
#argv <- add_argument(argv,"--bin",help="2, 8, 16 binned spots",default="16")
argv <- add_argument(argv,"--hvg",help="Number of highly variable genes to run PCA upon",default="2000")
argv <- add_argument(argv,"--resolution",help="clustering resolution",default=0.5)
argv <- add_argument(argv,"--min_umi",help="filter barcodes according to nCount",default=10)
#argv <- add_argument(argv,"--qs",help="The values of q to evaluate",default="2,10")
argv <- add_argument(argv, "--name", help = "the sample name", default = "spatial")
argv <- add_argument(argv, "--group", help = "the sample name", default = "spatial")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)

clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
clustcol <- c(color_protocol, clustcol)
col1 <- colorRampPalette(c("#7F0000","red","red","#FF7F00","#FF7F00","yellow","yellow","cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(clustcol)
corrcol <- colorRampPalette(c("red","orange","blue","white","white"))

#datapath = '/SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/test_HDspatial/20240412/Visium_HD_Human_Lung_Cancer/outs/'

#h5 <- "/SGRNJ06/randd/USER/zhouyiqi/work/analysis/space/mouse-brain-2/Mbrain_FFPE_96_96_5_1029/outs/filtered_feature_bc_matrix.h5"
#spatial <- '/SGRNJ06/randd/USER/zhouyiqi/work/analysis/space/mouse-brain-2/spatial/'
#outdir <- '/SGRNJ03/PiplineTest01/Pipline_test/fengshijing/spatial_dbit/mouse-brain-2/'
#name = 'mouse-brain-2'

name <- argv$name
group <- argv$group
h5 <- argv$h5
spatial <- argv$spatial

outdir <- argv$outdir
dir.create(outdir)
datapath = paste0(outdir,'/outs/')
dir.create(datapath)
system(paste('cp',h5,datapath))
system(paste0('cp -r ',spatial,' ',datapath,'/spatial/'))
hvg = as.numeric(argv$hvg)
pca = as.numeric(argv$pca)
resolution = as.numeric(argv$resolution)
min_umi <- as.numeric(argv$min_umi)


seu_data <- Load10X_Spatial(data.dir = datapath,image.scale = "hires")
seu_data$orig.ident <- name


#qs <- as.character(unlist(strsplit(argv$qs,split = ",")))
#name = 'human_lung'


assay <- 'Spatial'

meta <- seu_data@meta.data
ncount <- paste0('nCount_',assay)
cells_use <- rownames(meta[meta[[ncount]] > min_umi,])
print('Number of spots before filtering: ')
print(nrow(meta))

seu_data <- subset(seu_data, cells = cells_use)
print('Number of spots after filtering: ')
print(length(cells_use))

#prefix <- argv$prefix
# outdir = "/SGRNJ03/PiplineTest01/Pipline_test/fengshijing/Visium_HD/human_lung/"


####loading data####
#sce_data <- readVisium(datapath)
#seu_data <- Load10X_Spatial(data.dir = datapath, bin.size = bin)



###BayesSpace###
print("################### Star  Seurat Processing #############################")
###pre-processing data###
set.seed(519)

#bin_assay <- str_pad(as.character(bin),width = 3,pad = "0", side = 'left')

#if (bin == 2){assay = "Spatial.002um"}
#if (bin == 8){assay = "Spatial.008um"}
#if (bin == 16){assay = "Spatial.016um"}

print(assay)
DefaultAssay(seu_data) <- assay
seu_data$sample <- name
seu_data$group <- group
seu_data$sample <- factor(seu_data$sample)
seu_data$group <- factor(seu_data$group)

mito.genes <- grep(pattern = "^MT-", x=rownames(x=seu_data), value = TRUE, ignore.case = TRUE)
seu_data[["percent.mito"]] <- PercentageFeatureSet(seu_data, features = mito.genes)
Ribosomal <- grep(pattern = "^(RPL|RPS)", x = rownames(x=seu_data), value = TRUE, ignore.case = TRUE)
seu_data[["percent.Ribo"]] <- PercentageFeatureSet(seu_data, features = Ribosomal)


seu_data <- NormalizeData(seu_data)

seu_data <- FindVariableFeatures(seu_data, nfeatures = hvg)

seu_data <- ScaleData(seu_data)

seu_data <- RunPCA(seu_data, reduction.name = "pca")

seu_data <- FindNeighbors(seu_data, reduction = "pca", dims = 1:pca)

seu_data <- FindClusters(seu_data, resolution = resolution, cluster.name = "cluster")

seu_data <- RunUMAP(seu_data, reduction = "pca", reduction.name = "umap", dims = 1:pca)

seu_data$cluster <- factor(seu_data$cluster, levels = as.character(1:length(unique(seu_data$cluster))-1))
Idents(seu_data) <- 'cluster'
seu_data$raw_cluster <- seu_data$cluster


###1.QC###
qc_dir <- file.path(outdir, "1.QC")
dir.create(qc_dir)
  p1 <- VlnPlot(object = seu_data, features = c("percent.mito"), pt.size = 0,combine = TRUE, assay = assay, group.by = 'sample', cols = "greenyellow") + NoLegend() + labs(x = "")
  
#  p2 <- VlnPlot(object = seu_data, features = c("percent.Ribo"), pt.size = 0,combine = TRUE, assay = assay, cols = "darkorange",group.by = 'orig.ident') + NoLegend() + labs(x = "")
  
  p3 <- VlnPlot(object = seu_data, features = c(paste0("nFeature_",assay)), adjust = 1, pt.size = 0, combine = TRUE, assay = assay, cols = "darkslategray",group.by = 'sample') + labs(title = "nGene", x = "") + NoLegend()
  
  p4 <- VlnPlot(object = seu_data, features = c(paste0("nCount_",assay)), adjust = 1, pt.size = 0, combine = TRUE, assay = assay, cols = "orangered",group.by = 'sample') + labs(title = "nUMI", x = "") + NoLegend()
  pdf(paste0(qc_dir, '/', name, '_qc_vlnplot.pdf'), width = 8, height = 6)
  plot_grid(p1, p3, p4, ncol = 3)
  dev.off()
  png(paste0(qc_dir, '/', name, '_qc_vlnplot.png'), width = 8, height = 6,res=300,units='in')
  plot_grid(p1, p3, p4, ncol = 3)
  dev.off()


###2.Cluster
cluster_dir <- file.path(outdir, "2.Cluster")
dir.create(cluster_dir)

  pdf(paste0(cluster_dir, '/', name, '_PCElbowPlot.pdf'))
  ElbowPlot(seu_data, ndims = as.numeric(pca))
  dev.off()
  png(paste0(cluster_dir, '/', name, '_PCElbowPlot.png'))
  ElbowPlot(seu_data, ndims = as.numeric(pca))
  dev.off()




#####保存RDS文件
rds_dir <- file.path(outdir, "0.Rds")
dir.create(rds_dir)
saveRDS(seu_data, file = paste0(rds_dir, "/", name, ".rds"))
coords <- GetTissueCoordinates(seu_data)
write.table(coords,  paste0(rds_dir, "/spatial.xls"),sep = '\t',quote = F,row.names = T, col.names = NA)

#ident <- as.numeric(levels(x = seu_data))
#new <- ident + 1
#names(new) <- levels(x = seu_data)
#seu_data <- RenameIdents(seu_data, new)

clust <- summary(seu_data@active.ident)
cluster_spot <- as.data.frame(clust)
spot_number <- sum(cluster_spot$clust)
print(spot_number)
pt_use <- 0.6
if(spot_number > 1000){
  pt_use <- 0.4
}
if(spot_number > 2500){
  pt_use <- 0.3
}
if(spot_number > 4000){
  pt_use <- 0.2
}
if(spot_number > 5500){
  pt_use <- 0.15
}
if(spot_number > 6500){
  pt_use <- 0.1
}



  p1 <- FeaturePlot(seu_data, reduction = "umap", features = c("percent.mito"), pt.size = pt_use/2,raster = FALSE,order = T)+ scale_color_viridis()
  
  pdf(paste0(qc_dir, '/', name, '_percent.mito_umapFeatureplot.pdf'), width = 6, height = 6)
  print(p1)
  dev.off()
  png(paste0(qc_dir, '/', name, '_percent.mito_umapFeatureplot.png'), width = 6, height = 6,res=200,units='in')
  print(p1)
  dev.off()


  p3 <- FeaturePlot(seu_data, reduction = "umap", features =c(paste0("nFeature_",assay)), pt.size = pt_use/2, raster = FALSE) + labs(title = "nGene",order = T)+ scale_color_viridis()
  pdf(paste0(qc_dir, '/', name, '_nGene_umapFeatureplot.pdf'), width = 6, height = 6)
  print(p3)
  dev.off()
  png(paste0(qc_dir, '/', name, '_nGene_umapFeatureplot.png'), width = 6, height = 6,res=200,units='in')
  print(p3)
  dev.off()  
  
  

  p4 <- FeaturePlot(seu_data, reduction = "umap", features = c(paste0("nCount_",assay)),pt.size = pt_use/2, raster = FALSE, order = T) + labs(title = "nUMI")+ scale_color_viridis()
  pdf(paste0(qc_dir, '/', name, '_nUMI_umapFeatureplot.pdf'), width = 6, height = 6)
  print(p4)
  dev.off()
  png(paste0(qc_dir, '/', name, '_nUMI_umapFeatureplot.png'), width = 6, height = 6,res=200,units='in')
  print(p4)
  dev.off() 



  p1 <- SpatialFeaturePlot(seu_data, features = c("percent.mito"),  alpha = c(1, 1),crop = F, pt.size.factor = 1.5) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "percent.mito")
  pdf(paste0(qc_dir, '/', name, '_percent.mito_Tissueplot.pdf'), width = 6, height = 6)
  print(p1)
  dev.off()
  png(paste0(qc_dir, '/', name, '_percent.mito_Tissueplot.png'), width = 6, height = 6,res=300,units='in')
  print(p1)
  dev.off()  

  
  
  p3 <- SpatialFeaturePlot(seu_data, features = c(paste0("nFeature_",assay)),crop = F, pt.size.factor = 1.5, alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "nGene")
  pdf(paste0(qc_dir, '/', name, '_nGene_Tissueplot.pdf'), width = 6, height = 6)
  print(p3)
  dev.off()
  png(paste0(qc_dir, '/', name, '_nGene_Tissueplot.png'), width = 6, height = 6,res=300,units='in')
  print(p3)
  dev.off()    


  p4 <- SpatialFeaturePlot(seu_data, features = c(paste0("nCount_",assay)), crop = F, pt.size.factor = 1.5, alpha = c(1, 1)) + theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 15)) +labs(title = "nUMI")
  pdf(paste0(qc_dir, '/', name, '_nUMI_Tissueplot.pdf'), width = 6, height = 6)
  print(p4)
  dev.off()
  png(paste0(qc_dir, '/', name, '_nUMI_Tissueplot.png'), width = 6, height = 6,res=300,units='in')
  print(p4)
  dev.off()    
  


cluster_num <- data.frame(table(seu_data@active.ident))
colnames(cluster_num) <- c("cluster","cellnumber")
write.table(cluster_num,file=paste0(cluster_dir,'/',name,'_PercentPerSpot.xls'),sep='\t',quote=F,row.names=F)
label_value <- paste('(', round(cluster_num$cellnumber/sum(cluster_num$cellnumber) * 100, 1), '%)', sep = '')
label_plot <- paste(cluster_num$cluster, label_value, sep = '')
cluster_num$celltype <- label_plot
levels(cluster_num$cluster) <- levels(seu_data@active.ident)

p <- ggplot(data = cluster_num, mapping = aes(x = 'Content', y= cellnumber,fill = cluster)) + 
  geom_bar(stat = 'identity', position = 'stack',width=1,col='white') + coord_polar(theta = 'y') + 
  theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(panel.border=element_blank()) +
  labs(x = '', y = '', title = paste0('Percent of spots number per Cluster ', '(', name, ')')) + 
  scale_fill_manual(values=clustcol,name="celltype_prop",labels=sprintf("%s",cluster_num$celltype)) 

pdf(paste0(cluster_dir, '/', name, '_PercentPerSpot.pdf'))
p
dev.off()
png(paste0(cluster_dir, '/', name, '_PercentPerSpot.png'))
p
dev.off()

cluster.averages0 <- AggregateExpression(object = seu_data, return.seurat = TRUE, verbose = FALSE)
cluster.averages1 <- AggregateExpression(object = seu_data)
cluster.averages <- cluster.averages1[[assay]]
j <- length(levels(x = seu_data))
regions <- c()
for (i in 0: (j-1)){
  regions <- c(regions, paste0('cluster ', i))
}
colnames(cluster.averages) <- regions
write.table(cluster.averages, file= paste0(cluster_dir, '/', name, '_cluster_averages.xls'), sep = '\t', quote = FALSE, row.names = TRUE,col.names = NA)



#####降维图
spot_type <- as.data.frame(seu_data@active.ident)
colnames(spot_type) <- c("spot_type")
write.table(spot_type, file = paste0(cluster_dir, '/', name, '_spot_type.xls'), sep = '\t', quote = FALSE, row.names = TRUE,col.names = NA)

p <- DimPlot(seu_data, reduction = "pca",group.by='cluster', pt.size = pt_use, cols = clustcol,raster = F) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_pca.pdf'))
p
dev.off()
png(paste0(cluster_dir, '/', name, '_pca.png'))
p
dev.off()

p <- DimPlot(seu_data, reduction = "umap",group.by='cluster', pt.size = pt_use, cols = clustcol,raster = F) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_umap.pdf'))
p
dev.off()
png(paste0(cluster_dir, '/', name, '_umap.png'))
p
dev.off()

cluster_order = levels(seu_data@meta.data$cluster)
color_order = clustcol[1:length(cluster_order)]
cols <- list()
for(i in 1:length(cluster_order)) {
  cols[cluster_order[i]] <- color_order[i]
}
p <- SpatialDimPlot(seu_data,group.by='cluster', cols = cols,crop = F, pt.size.factor = 1.5) + labs(title = name) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(cluster_dir, '/', name, '_tissue.pdf'), width = 6, height = 6)
p
dev.off()
png(paste0(cluster_dir, '/', name, '_tissue.png'), width = 6, height = 6,res=300,units='in')
p
dev.off()


for (l in levels(seu_data)){

p1 <- SpatialDimPlot(seu_data,crop = F, pt.size.factor = 1.5, cells.highlight = WhichCells(seu_data, expression = cluster == l ), cols.highlight = c(cols[[l]], "lightgray")) + NoLegend() + ggtitle(paste("Cluster",l))

pdf(paste0(cluster_dir, '/', name,'.cluster',l, '_tissue.pdf'), width = 6, height = 6)
print(p1)
dev.off()
png(paste0(cluster_dir, '/', name,'.cluster',l, '_tissue.png'), width = 6, height = 6,res=300,units='in')
print(p1)
dev.off()

}


#####寻找差异基因####
diff_dir <- file.path(outdir, "3.Diff")
dir.create(diff_dir)
diffgene_dir <- file.path(diff_dir, "3.1.diffgene")
dir.create(diffgene_dir)

cluster_valid <- cluster_num[cluster_num$cellnumber > 3,]
markergenetop10 <- c()
for (l in cluster_valid$cluster){
  

  print(l)
  cluster.markers <- FindMarkers(object = seu_data, ident.1 = l, min.pct = 0.1, logfc.threshold = 0.25,max.cells.per.ident = 2000)
  colnames(cluster.markers)[2] <- c("avg_logFC")
  cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE),]
  write.table(data.frame(gene_id = rownames(cluster.markers),cluster.markers), file = paste0(diffgene_dir, '/', name, '_cluster', l, '_diffgenes.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
  top10 <- rownames(cluster.markers)[1:10]
  markergenetop10 <- c(markergenetop10,top10)
  print(top10)

}
seu_data.marker <- FindAllMarkers(object = seu_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
colnames(seu_data.marker)[2] <- c("avg_logFC")
markergene <- seu_data.marker %>% group_by(cluster) %>% top_n(2, avg_logFC)
markergenetop10 <- seu_data.marker %>% group_by(cluster) %>% top_n(10, avg_logFC)


diffheatmap_dir <- file.path(diff_dir, "3.2.heatmap")
dir.create(diffheatmap_dir)
heatname <- paste0(name, ' marker genes heatmap')


sub_seu <- subset(seu_data,downsample = 500)
sub_seu <- ScaleData(sub_seu, features =markergenetop10$gene)


pdf(paste0(diffheatmap_dir, '/', name, '_diffgene_DOHeatmapplot.pdf'))
DoHeatmap(sub_seu, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = FALSE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
dev.off()
png(paste0(diffheatmap_dir, '/', name, '_diffgene_DOHeatmapplot.png'))
DoHeatmap(sub_seu, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
dev.off()


heatname <- paste0(name, ' top10 marker genes heatmap')
pdf(paste0(diffheatmap_dir, '/', name, '_diffgenetop10_DOHeatmapplot.pdf'),width = 6,height = 8)
DoHeatmap(sub_seu, features = markergenetop10$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 6)) + scale_fill_viridis()
dev.off()
png(paste0(diffheatmap_dir, '/', name, '_diffgenetop10_DOHeatmapplot.png'),width = 6,height = 8,units = 'in', res = 200)
DoHeatmap(sub_seu, features = markergenetop10$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 6)) + scale_fill_viridis()
dev.off()

Allmarker <- seu_data.marker[order(seu_data.marker$avg_logFC, decreasing = TRUE),]
Allmarker <- Allmarker[,!colnames(Allmarker) %in% c("gene")]
write.table(Allmarker, file = paste0(diffgene_dir, '/', name, '_all_diffgenes.xls'), sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)

allmarker <- as.data.frame(markergene)
marker <- allmarker$gene
marker <- unique(marker)
print(marker)
pdf(paste0(diffheatmap_dir, '/', name, '_cluster_Heatmapplot.pdf'))
DoHeatmap(cluster.averages0, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = FALSE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
dev.off()
png(paste0(diffheatmap_dir, '/', name, '_cluster_Heatmapplot.png'))
DoHeatmap(cluster.averages0, features = markergene$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = FALSE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
dev.off()



######差异基因DotPlot
diffdotplot_dir <- file.path(diff_dir, "3.3.dotplot")
dir.create(diffdotplot_dir)

x_size <- 9
dot_size <- 7
if (length(marker) > 10){
  dot_size <- 6
}
if (length(marker) > 15){
  dot_size <- 6
}
if (length(marker) > 20){
  dot_size <- 5
}
if (length(marker) > 25){
  dot_size <- 4
}
if (length(marker) > 30){
  dot_size <- 3
}
x_size <- dot_size + 2.2
lab <- length(levels(seu_data))
l_size <- 10
if (lab > 10){
  l_size <- 9
}
if (lab > 15){
  l_size <- 8
}
if (lab > 20){
  l_size <- 7
}
if (lab > 25){
  l_size <- 6
}
if (lab > 30){
  l_size <- 5
}

pdf_h <- lab/2
pdf_w <- length(marker)/2

pdf(paste0(diffdotplot_dir, '/', name, '_diffgene_Dotplot.pdf'),width = pdf_w, height = pdf_h)
DotPlot(object = seu_data, features = marker, cols = c("blue", "red"), assay = assay, dot.scale = dot_size) + RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(color = guide_colorbar(title = "avg.exp.scale"), size = guide_legend(title = "pct.exp")) + theme(axis.text.x = element_text(size = x_size), axis.text.y = element_text(size = l_size))
dev.off()
png(paste0(diffdotplot_dir, '/', name, '_diffgene_Dotplot.png'),width = pdf_w, height = pdf_h,res = 300, units = 'in')
DotPlot(object = seu_data, features = marker, cols = c("blue", "red"), assay = assay, dot.scale = dot_size) + RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(color = guide_colorbar(title = "avg.exp.scale"), size = guide_legend(title = "pct.exp")) + theme(axis.text.x = element_text(size = x_size), axis.text.y = element_text(size = l_size))
dev.off()


###差异基因VlnPlot
diffvlnplot_dir <- file.path(diff_dir, "3.4.violinplot")
dir.create(diffvlnplot_dir)


for (l in as.character(unique(markergene$cluster))){

marker_l <- markergene[markergene$cluster == l,]

p1 <- VlnPlot(seu_data, features = marker_l$gene, pt.size = 0, cols = clustcol, ncol = 1, combine = TRUE) 
    pdf(paste0(diffvlnplot_dir, '/cluster', l, '_diffgene_vlnplot.pdf'))
    print(p1)
    dev.off()
    png(paste0(diffvlnplot_dir, '/cluster', l, '_diffgene_vlnplot.png'))
    print(p1)
    dev.off()


}


####差异基因FeaturePlot
difffeatureplot_dir <- file.path(diff_dir, "3.5.featureplot")
dir.create(difffeatureplot_dir)

for (l in as.character(unique(markergene$cluster))){

marker_l <- markergene[markergene$cluster == l,]

p1 <- FeaturePlot(seu_data, features = marker_l$gene,ncol = 2, combine = TRUE,   cols = c("lightgrey", "red"), pt.size = pt_use/2, raster = F) 
    pdf(paste0(difffeatureplot_dir, '/cluster', l, '_diffgene_umapFeatureplot.pdf'),width = 12,height = 6)
    print(p1)
    dev.off()
    png(paste0(difffeatureplot_dir, '/cluster', l, '_diffgene_umapFeatureplot.png'),width = 12,height = 6,res = 200, units = 'in')
    print(p1)
    dev.off()
p2 <- SpatialFeaturePlot(seu_data, marker_l$gene,ncol = 2, alpha = c(1, 1),crop = F, pt.size.factor = 1.5)
    pdf(paste0(difffeatureplot_dir, '/cluster', l, '_diffgene_Tissueplot.pdf'),width = 12,height = 6)
    print(p2)
    dev.off()
    png(paste0(difffeatureplot_dir, '/cluster', l, '_diffgene_Tissueplot.png'),width = 12,height = 6,res = 300, units = 'in')
    print(p2)
    dev.off()    

}




