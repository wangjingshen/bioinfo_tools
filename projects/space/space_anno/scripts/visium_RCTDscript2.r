#SeuratV5空转  sc&sp解卷

suppressMessages({
library(dplyr)
library(tidyr)
library(Seurat)
library(argparser)
library(gtools)
library(stringr)
library(spacexr)
library(arrow)
library(ggplot2)
library(viridis)
library(Cairo)
library(data.table)
library(presto)
# library(STdeconvolve)
# library(SPOTlight)
# library(NMF)
})

argv <- arg_parser('Pipeline of analysis of spatial transcriptomics with Seurat')
argv <- add_argument(argv, "--sc", help = "the path of the directory of single cell transcriptomics data", default = "")
argv <- add_argument(argv, "--sp", help = "the path of the directory of spatial transcriptomics data", default = "")
argv <- add_argument(argv, "--slot", help='Spatial.016um, Spatial.008um')
argv <- add_argument(argv, "--use_sample", help = "the sample name", default = "all")
# argv <- add_argument(argv, "--method", help = "", default = "rctd")
argv <- add_argument(argv, "--downsample", help = "", default = "all")
argv <- add_argument(argv, "--color_use", help = "the color of umap", default = "auto")
argv <- add_argument(argv, "--prefix", help = "", default = "prefix")
argv <- add_argument(argv, "--outdir", help = "", default = ".")
# argv <- add_argument(argv, "--mode", help='the doublet mode used for RCTD, "full" or "multi" mode is recommand for 10X Visiumi.', default='doublet')
# argv <- add_argument(argv, "--ct_method", help = "k.anchor", default = "bayes")
argv <- parse_args(argv)

############################################RCTD解卷
myRCTD <- function(scrds,sprds){
	##scrds处理——reference
    # scrds@meta.data$cluster <- scrds@active.ident
	# scrds@meta.data <- scrds@meta.data %>% mutate(cluster = str_replace(cluster,"/","_"))
	# Idents(scrds) <- 'cluster'
	a <- as.data.frame(table(Idents(scrds))) %>% filter(Freq < 25) #%>% select(-Freq)
	if(nrow(a)>0){scrds <- subset(x = scrds, idents = a[['Var1']], invert = TRUE)}
	
	counts <- scrds@assays$RNA@counts
	cell_types <- select(scrds@meta.data, cluster)
	cell_types$barcodes <- rownames(cell_types)
	cell_types <- setNames(cell_types[[1]], cell_types[[2]]) #
	cell_types <- as.factor(cell_types) 
	nUMI <- select(scrds@meta.data, nCount_RNA)
	nUMI$barcodes <-  rownames(nUMI)
	nUMI <- setNames(nUMI[[1]], nUMI[[2]])
	reference <- Reference(counts, cell_types, nUMI,require_int = FALSE)
	print(dim(reference@counts))
	table(reference@cell_types)
	##sprds处理——puck
	counts <- sprds@assays[[slot]]@layers$counts
	rownames(counts) <- rownames(sprds)
	colnames(counts) <- colnames(sprds)
	
	coords <- GetTissueCoordinates(sprds)
	colnames(coords) <- c("x", "y")
	nUMI <- colSums(counts)
	names(nUMI) <- colnames(sprds)
	puck <- SpatialRNA(coords, counts, nUMI)
	myRCTD <- create.RCTD(puck, reference, max_cores = 2)
	myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
	return(myRCTD)
}


##绘图
sp_plot <- function(sprds,outdir,prefix,clustcol){

	p <- SpatialDimPlot(sprds, cols = clustcol[levels(sprds)], group.by = "celltype")
	# ggsave(paste0(outdir,'/',prefix,'_cell_type.pdf'))
	# ggsave(paste0(outdir,'/',prefix,'_cell_type.png'))
	CairoPDF(paste0(outdir,'/',prefix,'_cell_type.pdf'))
	print(p)
	dev.off()
	
	CairoPNG(paste0(outdir,'/',prefix,'_cell_type.png'))
	print(p)
	dev.off()
	

	p_pca <- DimPlot(sprds,cols = clustcol[levels(sprds)],reduction = 'pca')
	ggsave(paste0(outdir,'/',prefix,".pca.pdf"))
	ggsave(paste0(outdir,'/',prefix,".pca.png"))

	p_dim <- DimPlot(sprds,cols = clustcol[levels(sprds)],reduction = 'umap')
	ggsave(paste0(outdir,'/',prefix,".umap.pdf"))
	ggsave(paste0(outdir,'/',prefix,".umap.png"))
	
	for(ct in levels(sprds)){
		sub_sprds <- subset(sprds, celltype == ct)
		if(dim(sub_sprds)[2] < 2) next

		fea <- SpatialFeaturePlot(sub_sprds, features='weight')
		# ggsave(paste0(outdir,'/weights/',ct,"_weight.pdf"))
		# ggsave(paste0(outdir,'/weights/',ct,"_weight.png"))
		CairoPDF(paste0(outdir,'/weights/',ct,"_weight.pdf"))
		print(fea)
		dev.off()
	
		CairoPNG(paste0(outdir,'/weights/',ct,"_weight.png"))
		print(fea)
		dev.off()
	}

}

cellprop_plot <- function(sprds,outdir,prefix,clustcol){

	meta <- sprds@meta.data

	tbl <- table(meta$seurat_clusters,meta$celltype)
	tbl2 <- apply(tbl,1,function(x)x/sum(x))
	tbl3 <- as.data.frame(t(tbl2))
	data.m <- tbl3 %>% pivot_longer(everything(),names_to='Celltype',values_to = "Fraction") %>% as.data.frame()
	rowname <- rep(rownames(tbl3),each=length(levels(sprds)))
	data.m$Cluster <- rowname

	q <- length(unique(data.m$Cluster))-1
	data.m$Cluster <- factor(data.m$Cluster, levels = rev(0:q))
	
	p <- ggplot(data.m,aes(x=Cluster,y=Fraction,fill=Celltype))+
	geom_bar(stat="identity",width=0.7,color = "black")+
	theme_classic()
	
	p <- p + scale_fill_manual(values = clustcol[levels(sprds)]) +
	xlab ('') + ylab ('Fraction of Cells') + ylim(c(0,1))+
	scale_y_continuous(expand=c(0,0))+
	theme (axis.line.y = element_line(colour = "transparent"),
		   axis.ticks.y = element_blank(),
		   axis.text.y = element_text(angle = 0,color = "black"),
		   aspect.ratio=1) +
	guides(fill=guide_legend(title='',ncol = 1)) + coord_flip()

	write.table(tbl3, paste0(outdir,'/',prefix,'.CelltypePerCluster.xls'), sep = '\t',quote = F, row.names = T, col.names = NA)

	pdf(paste0(outdir,'/',prefix,'.CelltypePerCluster.pdf'),width = 7,height = 7)
	print(p)
	dev.off()
	png(paste0(outdir,'/',prefix,'.CelltypePerCluster.png'),width = 7,height = 7,res = 300,units = 'in')
	print(p)
	dev.off()
}

# 差异基因
diffgene <- function(sprds,assay,outdir,prefix,clustcol){

	diff_dir <- file.path(outdir, "diffgene")
	dir.create(diff_dir)

	cluster_num <- data.frame(table(sprds@active.ident))
	colnames(cluster_num) <- c("celltype","cellnumber")
	cluster_valid <- cluster_num[cluster_num$cellnumber > 3,]

	for (ct in cluster_valid$celltype){
	  print(ct)
	  cluster.markers <- FindMarkers(object = sprds, ident.1 = ct, min.pct = 0.1, logfc.threshold = 0.25,max.cells.per.ident = 2000)
	  colnames(cluster.markers)[2] <- c("avg_logFC")
	  cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE),]
	  write.table(data.frame(gene_id = rownames(cluster.markers),cluster.markers), file = paste0(diff_dir, '/', prefix, '_cluster', ct, '_diffgenes.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
	  top10 <- rownames(cluster.markers)[1:10]
	  print(top10)
	}
	seu_data.marker <- FindAllMarkers(object = sprds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	colnames(seu_data.marker)[2] <- c("avg_logFC")

	markergenetop10 <- seu_data.marker %>% group_by(cluster) %>% top_n(10, avg_logFC)

	diffheatmap_dir <- file.path(diff_dir, "heatmap")
	dir.create(diffheatmap_dir)

	sub_seu <- subset(sprds,downsample = 500)
	sub_seu <- ScaleData(sub_seu, features =markergenetop10$gene)

	pdf(paste0(diffheatmap_dir, '/', prefix, '_diffgenetop10_DOHeatmapplot.pdf'),width = 9,height = 16)
	print(DoHeatmap(sub_seu, features = markergenetop10$gene, group.colors=clustcol[levels(sub_seu)], size = 3, vjust = 1, angle = 90, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 6)) + scale_fill_viridis())
	dev.off()
	
	png(paste0(diffheatmap_dir, '/', prefix, '_diffgenetop10_DOHeatmapplot.png'),width = 9,height = 16,units = 'in', res = 200)
	print(DoHeatmap(sub_seu, features = markergenetop10$gene, group.colors=clustcol[levels(sub_seu)], size = 3, vjust = 1, angle = 90, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = assay) + theme(axis.text.y = element_text(size = 6)) + scale_fill_viridis())
	dev.off()
}

############主流程
##读取数据
outdir <- argv$outdir
prefix <- argv$prefix
# ct_method <-argv$ct_method
dir.create(outdir, recursive = TRUE)
resdir <- paste0(outdir,'/',prefix)
dir.create(resdir)
dir.create(paste0(resdir,'/weights'))
sprds <- readRDS(argv$sp)
slot <- argv$slot
scrds <- readRDS(argv$sc)
# scrds <- UpdateSeuratObject(scrds)
scrds@meta.data$cluster <- scrds@active.ident
scrds@meta.data <- scrds@meta.data %>% mutate(cluster = str_replace(cluster,"/","_"))
Idents(scrds) <- 'cluster'

a <- as.data.frame(table(Idents(scrds))) %>% filter(Freq < 25) #%>% select(-Freq)
if(nrow(a)==0){scrds2 <- scrds}else{scrds2 <- subset(x = scrds, idents = a[['Var1']], invert = TRUE)}

if(argv$downsample != 'all'){
	ncell = as.numeric(argv$downsample)
	scrds <- subset(scrds, downsample=ncell)
}else{
	scrds = scrds
}

if(argv$use_sample != 'all'){
	use_sample = unlist(strsplit(argv$use_sample,split=','))
	scrds <- subset(scrds,cells=rownames(scrds@meta.data)[which(scrds$sample %in% use_sample)])
}

if('cluster_colors' %in% colnames(scrds@meta.data)){
	out_clustcol <- data.frame(scrds@meta.data$cluster,scrds@meta.data$cluster_colors)
	out_clustcol <- unique(out_clustcol)
	colnames(out_clustcol) <- c('cluster','cmap')
	clustcol <- out_clustcol$cmap
	cts <- out_clustcol$cluster
	names(clustcol) <- cts
}else{
	if(argv$color_use != 'auto'){
		clustcol <- argv$color_use
	}else{
	  color_use <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#C6CCC3","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101","#9370DB","#EE9A00","#27408B")
	  source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
	  clustcol <- c(color_protocol, color_use)
	}
	cts <- levels(scrds)
	clustcol <- clustcol[1:length(cts)]
	names(clustcol) <- cts
	# out_clustcol <- data.frame(clustcol[1:length(cts)])
	# out_clustcol <- cbind(rownames(out_clustcol),out_clustcol)
	# colnames(out_clustcol) <- c('cluster','cmap')
}

###解卷积
print("start deconvolution")

rctd <- myRCTD(scrds,sprds)
saveRDS(rctd,paste0(outdir,'/',prefix,".RCTD.rds"))

weight <- rctd@results$weights_doublet
colnames(weight) <- c('first_type_weights','second_type_weights')
prediction <- cbind(rctd@results$results_df[,c('spot_class','first_type','second_type')], weight)

prediction <- prediction %>% mutate(celltype = as.character(first_type), celltype_colors= clustcol[as.character(first_type)],weight = first_type_weights)

print("deconvolution over")

sprds <- subset(sprds, cells = rownames(prediction))
sprds <- AddMetaData(sprds, metadata = prediction)
# sprds@meta.data$celltype_prediction <- gsub('-','_',sprds@meta.data$celltype_prediction)
Idents(sprds) <- sprds$celltype

ct_level <- as.character()
for(i in 1:length(levels(scrds))){
	if(levels(scrds)[i] %in% levels(sprds))
		ct_level <- append(ct_level,levels(scrds)[i])
}
levels(sprds) <- ct_level

saveRDS(sprds,paste0(outdir,'/',prefix,".spatial.rds"))

p <- SpatialDimPlot(sprds, group.by = "spot_class")
ggsave(paste0(outdir,'/',prefix,'_spot_class.pdf'))
ggsave(paste0(outdir,'/',prefix,'_spot_class.png'))

sp_plot(sprds,resdir,prefix,clustcol)
cellprop_plot(sprds,resdir,prefix,clustcol)

diffgene(sprds,slot,resdir,prefix,clustcol)

