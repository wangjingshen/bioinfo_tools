suppressMessages(suppressWarnings({
    library(Seurat)
    library(harmony)
    library(argparser)
    library(tidyverse)
    library(patchwork)
    library(future)
    library(furrr)
    library(SingleR)
    library(logger)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--matrix_10X", help = "10X matrix, split by ,")
argv <- add_argument(argv, "--spname", help = "sample name, split by ,")
argv <- add_argument(argv, "--gname", help = "group name, split by ,")
argv <- add_argument(argv, "--rm_batch", help="rm batch, T or F. Default: F")
argv <- add_argument(argv, "--rm_batch_var", help="rm batch var. Default: sample")
argv <- add_argument(argv, "--resolution", help="resolution. Default: 0.8")
argv <- add_argument(argv, "--species", help="species, for annotation. Default: human")
argv <- add_argument(argv, "--gtf_biotype", help="gtf_biotype tsv")
argv <- add_argument(argv, "--gene_type", default = "totalRNA", help="gene_type. Default: totalRNA")
argv <- add_argument(argv, "--outdir", default = "outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

matrix_10X <- str_split(argv$matrix_10X, ",", simplify = TRUE)
spname <- str_replace_all(str_split(argv$spname, ",", simplify = TRUE), "-", "_")
#spname <- gsub("-","_",unlist(strsplit(argv$spname, split = ",")))
gname <- str_replace_all(str_split(argv$gname, ",", simplify = TRUE), "-", "_")

if (length(unique(c(length(matrix_10X), length(spname), length(gname)))) > 1) {
  stop("The quantities of matrix_10X, spname and gname must be consistent！")
}

rm_batch <- ifelse(is.na(argv$rm_batch), "F", argv$rm_batch)
rm_batch_var <- ifelse(is.na(argv$rm_batch_var), "sample", argv$rm_batch_var)
resolution <- ifelse(is.na(argv$resolution), 0.8, as.numeric(argv$resolution))
species <- ifelse(is.na(argv$species), "human", argv$species)
gtf_biotype <- argv$gtf_biotype
gene_type <- ifelse(is.na(argv$gene_type), "totalRNA", argv$gene_type)

outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
plot_dir <- str_glue("{outdir}/plot/")
if(!dir.exists(plot_dir)){
    dir.create(plot_dir, recursive = TRUE)
}

# read data --
data_seurat_list <- pmap(list(matrix_10X, spname, gname),
    function(dir, sp, gp){
        data <- Read10X(dir) 
        if(gene_type == "protein_coding"){
            cat("pretein_coding subset...\n")
            gtf_biotype = read.table(gtf_biotype, sep="\t")
            gene_sub <- intersect(row.names(data), gtf_biotype$V2[ gtf_biotype$V3 == "protein_coding"])
            data <- data[gene_sub, ]
            cat("data shape: \n")
            print(dim(data))
        }
        if(gene_type == "non_coding"){
            cat("non_coding subset...\n")
            gtf_biotype = read.table(gtf_biotype, sep="\t")
            gene_sub <- intersect(row.names(data), gtf_biotype$V2[ gtf_biotype$V3 != "protein_coding"])
            data <- data[gene_sub, ]
            cat("data shape: \n")
            print(dim(data))         
        }
        
        data <- CreateSeuratObject(counts = data)
        data@meta.data$sample <- sp
        data@meta.data$group <- gp
        return(data)
    }
)
names(data_seurat_list) <- spname
# merge data --
if(length(spname)==1){
    data_seurat <- data_seurat_list[[1]]
    data_seurat <- RenameCells(data_seurat, add.cell.id = names(data_seurat_list))
}else{
    data_seurat <- merge(data_seurat_list[[1]], y = data_seurat_list[-1], add.cell.ids = names(data_seurat_list))
}

if(species %in% c("human", "mouse")){
    data_seurat$percent.mt <- ifelse(species == "human", 
                                     PercentageFeatureSet(data_seurat, pattern = "^MT-"),
                                     PercentageFeatureSet(data_seurat, pattern = "^Mt-"))
}

data_seurat <- data_seurat %>% 
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData()%>% 
    RunPCA(verbose =F) 

# clustering --
if(rm_batch %in% c("T", "True", "TRUE")){
    my_harmony_embeddings <- HarmonyMatrix(
        data_mat  = as.matrix(Embeddings(data_seurat)),
        meta_data = data_seurat@meta.data,
        vars_use  = rm_batch_var,
        do_pca = FALSE)

    rownames(my_harmony_embeddings) <- rownames(Embeddings(data_seurat))
    data_seurat[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(data_seurat))
}

data_seurat <- data_seurat %>%
    FindNeighbors(reduction = ifelse(rm_batch %in% c("T", "True", "TRUE"), "harmony", "pca"), dims = 1:20, verbose = F) %>%
    FindClusters(resolution = resolution, verbose = F) %>%
    RunTSNE(reduction = ifelse(rm_batch %in% c("T", "True", "TRUE"), "harmony", "pca"), dims = 1:20, do.fast = TRUE, check_duplicates = FALSE, verbose = F) %>%
    RunUMAP(reduction = ifelse(rm_batch %in% c("T", "True", "TRUE"), "harmony", "pca"), dims = 1:20, verbose = F)

# plot
function_save_plot <- function(p, name, h, w){
    ggsave(str_glue("{plot_dir}/{name}.png"), plot = p, height = h, width = w)
    ggsave(str_glue("{plot_dir}/{name}.pdf"), plot = p, height = h, width = w) 
}

dim_plots <- list(
    sample = DimPlot(data_seurat, group.by = "sample", label = T, repel = T),
    group = DimPlot(data_seurat, group.by = "group", label = T, repel = T),
    seurat_clusters = DimPlot(data_seurat, group.by = "seurat_clusters", label = T, repel = T)
)
future_walk2(names(dim_plots), dim_plots,
            function(name, p){
                function_save_plot(p = p, name = name, h = 6.5, w = 8)
            }, .options = furrr_options(seed = NULL))

# singleR
if(species %in% c("human", "mouse")){
    ref_path <- switch(
        species,
        human = "/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/HumanPrimaryCellAtlasData.rds",
        mouse = "/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/MouseRNAseqData.rds"
    )
    ref = readRDS(ref_path)
    anno <- SingleR(test = GetAssayData(data_seurat, slot = "data"), ref = ref, 
                    clusters = data_seurat$seurat_clusters,
                    assay.type.test=1, labels = ref$label.main)
    data_seurat$cluster <- plyr::mapvalues(x = data_seurat$seurat_clusters, from = row.names(anno), to = anno$labels) %>%
        factor(levels = sort(unique(.)))
    function_save_plot(DimPlot(data_seurat, group.by = "cluster", label = TRUE), "seurat_cluster_singleR", h = 6.5, w = 8)
}

# cellular_composition
p2 <- ggplot(data_seurat@meta.data, aes(x=sample, fill=cluster))+
    geom_bar(color="black",position="fill",width = 0.5)+
    theme_bw()+
    theme(legend.title = element_blank(), axis.title = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))
ggsave(str_glue("{plot_dir}/cellular_composition.png"), plot = p2, height = 6, width = 5)
ggsave(str_glue("{plot_dir}/cellular_composition.pdf"), plot = p2, height = 6, width = 5)

write.table(
    xtabs(~ sample + cluster, data_seurat@meta.data) %>% as.data.frame.matrix(),
    str_glue("{outdir}/sample_cluster.xls"), quote = F, sep = "\t", col.names = F)
#write.table(t(as.data.frame(table(data_seurat$sample,data_seurat$cluster)) %>% spread(Var2, Freq) ),
#         str_glue("{outdir}/sample_cluster.xls"),quote=F,sep="\t",col.names=F)

# saveRDS --
saveRDS(data_seurat, str_glue("{outdir}/seurat.rds"))
