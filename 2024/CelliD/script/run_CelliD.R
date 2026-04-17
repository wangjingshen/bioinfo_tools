## CelliD 

# load packages --
suppressWarnings(suppressMessages({
  library(argparser)
  library(Seurat)
  library(CellID)
  library(tidyverse) 
  library(ggpubr)
}))

# args ----
args <- arg_parser('')
args <- add_argument(args,"--rds", help="seurat rds file")
args <- add_argument(args,"--resolution", help="Optional, resolution, Default: 'seurat_clusters' in rds")
args <- add_argument(args,"--species", help="Optional,species, Hs(Default, human) or Mm(mouse)")
args <- add_argument(args,"--mode", help="analysis mode, all(markers of all organ) or single(markers of single organ) or common(markers of single organ and common celltype)")
args <- add_argument(args,"--organ", help="name of organ, raw reference: Adrenal glands, Blood, Bone, Brain, Connective tissue, Embryo, Epithelium, Eye, GI tract, Heart, Immune system, Kidney, Liver,
                     Lungs, Mammary gland, Olfactory system, Oral cavity, Pancreas, Placenta, Reproductive, Skeletal muscle, Skin, Smooth muscle, Thymus, Thyroid, Urinary bladder, Vasculature, Zygote,
                     update reference: BM_PBMC")
args <- add_argument(args,"--ref",help="Optional, markers reference, raw(Default, PanglaoDB_markers) or update(markers from data scienece)")
args <- add_argument(args,"--nfeatures",help="Optional, integer of top n features to consider for hypergeometric test, Default: 200, larger number(for example, 700) may be suitable for update ref")
args <- add_argument(args,"--outdir",help="Optional, outdir, Default: outdir")
args <- parse_args(args)

resolution <- ifelse(is.na(args$resolution), "seurat_clusters", paste0("RNA_snn_res.",args$resolution))
species <- ifelse(is.na(args$species), "Hs", args$species)
nfeatures <- ifelse(is.na(args$nfeatures), 200, as.numeric(args$nfeatures))
reference <- ifelse(is.na(args$ref), "raw", args$ref)
outdir <- ifelse(is.na(args$outdir), "outdir", args$outdir)

# marker reference --
if(reference == "raw"){
  ref <- read_tsv("/SGRNJ06/randd/USER/wangjingshen/project/CelliD/data/PanglaoDB_markers_27_Mar_2020.tsv")
}
if(reference == "update"){
  ref <- read_tsv("/SGRNJ06/randd/USER/wangjingshen/project/CelliD/data/markers_update.tsv")
}

# create outdir
if(!file.exists(outdir)){
  dir.create(outdir)
}


# read data and make ref gene set ----
data_seurat <- readRDS(args$rds)
data_seurat <- RunMCA(data_seurat)

# make ref gene set --
common_cells <- c("B cells","Basophils","Dendritic cells","Eosinophils","Macrophages","Mast cells",
                  "Monocytes","Natural killer T cells","Neutrophils","NK cells","Plasma cells",
                  "Plasmacytoid dendritic cells","T cells",
                  "Basal cells","Endothelial cells","Epithelial cells","Fibroblasts","Myofibroblasts",
                  "Pericytes","Smooth muscle cells","Stromal cells")

if(args$mode == "all"){
  ref_subset <- ref %>%  
    filter(str_detect(species, species))%>%  
    group_by(`cell type`) %>%  
    summarise(geneset = list(`official gene symbol`))
  ref_gs <- setNames(ref_subset$geneset, ref_subset$`cell type`)
}
if(args$mode == "single"){
  ref_subset <- ref %>% 
    filter(organ %in% unlist(str_split(args$organ,pattern = ",",simplify = F))) %>%  
    filter(str_detect(species, species)) %>%  
    group_by(`cell type`) %>%  
    summarise(geneset = list(`official gene symbol`))
  ref_gs <- setNames(ref_subset$geneset, ref_subset$`cell type`)
}
if(args$mode == "common"){
  ref_subset <- ref %>% 
    filter(organ %in% unlist(str_split(args$organ,pattern = ",",simplify = F))  | `cell type` %in% common_cells) %>%  
    filter(str_detect(species, species)) %>%  
    group_by(`cell type`) %>%  
    summarise(geneset = list(`official gene symbol`))
  ref_gs <- setNames(ref_subset$geneset, ref_subset$`cell type`)
}

# prediction ----
HGT <- RunCellHGT(data_seurat, pathways = ref_gs, dims = 1:50, n.features = nfeatures, minSize = 1)
prediction <- rownames(HGT)[apply(HGT, 2, which.max)]   # max -log10 corrected p-value
data_seurat$CelliD <- prediction
# prediction_signif <- ifelse(apply(HGT, 2, max)>2, yes = prediction, "unassigned")   # -log10p >2 
#data_seurat$CelliD_signif <- prediction_signif

# annotate cluster --
data_seurat$cluster_CelliD <- plyr::mapvalues(data_seurat@meta.data[[resolution]], 
    from = unique(data_seurat@meta.data[[resolution]]), 
    to = sapply(unique(data_seurat@meta.data[[resolution]]), function(x){
                names(which.max(table(data_seurat$CelliD[ data_seurat@meta.data[[resolution]] == x])))}))
p <- DimPlot(data_seurat, group.by = "cluster_CelliD", label = T, label.size = 5, repel = T)
ggsave(paste0(outdir, "/seurat_cluster_CelliD.pdf"), plot=p, height = 6.5, width = 8)
ggsave(paste0(outdir, "/seurat_cluster_CelliD.png"), plot=p, height = 6.5, width = 8)


# cellular_composition
p2 <- data_seurat@meta.data %>%
    ggplot(aes(x=sample,fill=cluster_CelliD))+
        geom_bar(color="black",position="fill",width = 0.5)+
        theme_bw()+
        theme(legend.title = element_blank(),axis.title = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
ggsave(paste0(outdir, "/cellular_composition.png"), plot= p2, width = 4, height = 5)
ggsave(paste0(outdir, "/cellular_composition.pdf"), plot= p2, width = 4, height = 5)


function_freq_stat <- function(cluster_id){
    freq <- as.data.frame(table(data_seurat$CelliD[data_seurat@meta.data[[resolution]] == cluster_id])/sum(data_seurat@meta.data[[resolution]] == cluster_id))
    freq <- freq[order(freq$Freq, decreasing = T),]
    freq <- freq[ freq$Freq > 0.05,]
    freq$Var1 <- as.character(freq$Var1)
    freq[nrow(freq)+1, 1]="others"
    freq[nrow(freq), 2] = 1- sum(freq[1:nrow(freq)-1, 2])
    freq["cluster"] = cluster_id
    colnames(freq) <- c("cell_type", "freq", "cluster")
    return(freq[, c(3,1,2)])
}
freq_stat <- do.call(rbind,lapply(sort(as.numeric(as.character(unique(data_seurat@meta.data[[resolution]])))), 
                                       function_freq_stat))
#write.table(freq_stat, paste0(outdir, "/CelliD_prediction_freq_each_cluster.tsv"), sep="\t", quote=F, row.names = F)

# output CelliD prediction --
#write.table(data.frame( row.names = colnames(data_seurat),
#                          cell = colnames(data_seurat),
#                          cluster_CelliD = data_seurat$cluster_CelliD), 
#paste0(outdir, "/cluster_CelliD.txt"), sep="\t", quote=F, row.names = F)


saveRDS(data_seurat, paste0(outdir, "/seurat_CelliD.rds"))


# end ----
cat("#-------------\n")
cat("CelliD done.")