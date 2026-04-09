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
argv <- add_argument(argv, "--rds", help = "Seurat rds")
argv <- add_argument(argv, "--spname", help = "spname")
argv <- add_argument(argv, "--anno_var", help = "anno_var")
argv <- add_argument(argv, "--marker_ref", help = "marker reference tsv")
argv <- parse_args(argv)

spname <- argv$spname
rds <- argv$rds
anno_var <- ifelse(is.na(argv$anno_var), "seurat_clusters", argv$anno_var)
marker_ref <- argv$marker_ref

#create dir
out_dir = str_glue('{spname}/auto_assign/')
plot_dir = str_glue('{spname}/auto_assign/plot/')
if(!dir.exists(plot_dir)){
    dir.create(plot_dir, recursive = TRUE)
}

# read rds
rds <- readRDS(rds)
marker_file <- read_tsv(marker_ref)
ref_cluster <- (marker_file)[, 1, drop = T]

seurat_clusters <- sort(unique(rds$seurat_clusters))
#seurat_clusters <- sort(unique(rds@active.ident))

# type marker
flag = 0 
for (cluster in seurat_clusters){
	index = 0
	for (celltype in ref_cluster){
    	index = index + 1
    	pos = unlist(strsplit(marker_file[index, 2, drop = T], ","))
    	neg = tryCatch(unlist(strsplit(marker_file[index, 3, drop=T], ",")), error=function(e){} )
    	for (F in pos){
      		tryCatch({
        		diff_sub <- FindMarkers(rds, feature=F, ident.1=cluster, min.pct = 0, logfc.threshold = -Inf)
        		diff_sub$cell_type <- celltype
        		diff_sub$cluster <- cluster
        		diff_sub <- rownames_to_column(diff_sub, var="gene")
        		diff_sub$type <- "positive"
        		if (flag==0){
          			df_diff <- diff_sub
          			flag = flag + 1
        		} else {
          			df_diff <- rbind(df_diff, diff_sub)
          		}
        	}
        	,error=function(e){print(paste0(F," not found in cluster ",cluster)) })
    	}

    	if (!is.na(neg) && !is.null(neg)){
    		for (F in neg){
      			tryCatch({
        			diff_sub <- FindMarkers(rds, feature=F, ident.1=cluster, min.pct = 0, logfc.threshold = -Inf)
        			diff_sub$cell_type <- celltype
        			diff_sub$cluster <- cluster
        			diff_sub <- rownames_to_column(diff_sub,var="gene")
        			diff_sub$type <- "negative"
        			if (flag==0){
          				df_diff <- diff_sub
          				flag = flag + 1
        			} else {
          				df_diff <- rbind(df_diff, diff_sub)
          			}
        		}
        	,error=function(e){print(paste0(F," not found in cluster ",cluster)) })
    		}
    	}
	}
}

#df_diff$cluster <- as.numeric(df_diff$cluster) + 1    # old version
df_diff <- mutate(df_diff, pct.diff = pct.1 - pct.2)
write_tsv(df_diff, str_glue('{out_dir}/{spname}_type_marker_exp.tsv'))


# plot
df_diff_gb <- group_by(df_diff, cluster, cell_type)
for (cluster in seurat_clusters){
	df_diff_gb_sub = df_diff_gb[df_diff_gb$cluster==cluster,]

	p1 <- ggplot(df_diff_gb_sub, aes(x = interaction(gene, cell_type, type), y=pct.diff, fill=cell_type)) +
    	geom_bar(stat = "identity") +
    	coord_flip() +
    	scale_fill_manual(values = color_protocol)
	ggsave(str_glue('{plot_dir}/{cluster}_pctdiff.png'), plot = p1, width=12, height=10)

	p2 <- ggplot(df_diff_gb_sub, aes(x = interaction(gene, cell_type, type), avg_log2FC, fill=cell_type)) +
    	geom_bar(stat = "identity") +
    	coord_flip() +
    	scale_fill_manual(values = color_protocol)
	ggsave(str_glue('{plot_dir}/{cluster}_logfc.png'), plot = p2, width=12, height=10)
}

# auto assign
df_diff[df_diff$type=="negative",]$avg_log2FC = -(df_diff[df_diff$type=="negative",]$avg_log2FC)
df_diff[df_diff$type=="negative",]$pct.diff = -(df_diff[df_diff$type=="negative",]$pct.diff)
df_diff_gb <- df_diff %>%
	group_by(cluster, cell_type) %>%
	summarize(avg_pct.diff = mean(pct.diff),
			  avg_log2FC = mean(avg_log2FC),
			  max_p_val_adj = max(p_val_adj))

df_anno <- df_diff_gb %>%
	ungroup() %>% 
	group_by(cluster) %>%
	mutate(pct_rank = rank(avg_pct.diff),
           logfc_rank = rank(avg_log2FC),
		   total_rank = pct_rank + logfc_rank) %>%
	filter(total_rank == max(total_rank)) %>%
	arrange(as.numeric(cluster)) %>%
	select(cluster, cell_type, avg_pct.diff, avg_log2FC, max_p_val_adj)

#as2 <- as1 %>% ungroup %>% 
#	group_by(cluster) %>% 
#	filter(total_rank == max(total_rank)) %>%
#	arrange(as.numeric(cluster))

#as3 <- select(as2, cluster, cell_type, avg_pct.diff, avg_log2FC, max_p_val_adj)
df_anno[( df_anno$avg_pct.diff < 0 | df_anno$avg_log2FC < 0),]$cell_type = 'NA'
write_tsv(df_anno, stringr::str_glue('{out_dir}/{spname}_auto_cluster_type.tsv'))