## code to visualize results of scFEA ----

cat("\n*** 03.scFEA visualization start running...\n")
timestart <- proc.time()
  
# library packages --
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(pheatmap))
suppressMessages(library(argparser))

# args --
args <- arg_parser('')
args <- add_argument(args,"--rds", help="seurat rds file")
args <- add_argument(args,"--cell_type", help="name of variable about cell type")
args <- add_argument(args,"--group", help="name of variable about group")
args <- add_argument(args,"--model", help="model, all or split")
args <- parse_args(args)

rds <- readRDS(args$rds)
model <- args$model
cell_type <- args$cell_type
group <- args$group
module_info <- read.table("/SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/data/module_info_add_supermodule.tsv",sep="\t")
 
# barplot  --
function_barplot <- function(model, flux, balance, metadata_sub, cluster){
    cluster <- gsub(" ", "-", cluster)
    metadata_sub_group <- metadata_sub[[group]]
    metadata_sub_cell_type <- metadata_sub[[cell_type]]
    
    # set group and name --
    if(model == "split"){
        g1 <- row.names(metadata_sub)[metadata_sub_group == unique(metadata_sub_group)[1]]
        g2 <- row.names(metadata_sub)[metadata_sub_group == unique(metadata_sub_group)[2]] 
        name = ""
    }
    if(model == "all"){
        g1 <- row.names(metadata_sub)[metadata_sub_group == unique(metadata_sub_group)[1] & metadata_sub_cell_type == cluster]
        g2 <- row.names(metadata_sub)[metadata_sub_group == unique(metadata_sub_group)[2] & metadata_sub_cell_type == cluster] 
        name = "all_"
    }
    ## flux --
    flux_barplot <- data.frame(module_name = factor(row.names(flux),levels = row.names(flux)),
                               flux_g1 = apply(flux[,g1], 1, mean),
                               flux_g2 = apply(flux[,g2], 1, mean),
                               supermodule_name = word(row.names(flux),1,sep=fixed(":")),
                               flux_change = apply(flux[,g1], 1, mean) / apply(flux[,g2], 1, mean)) 
    # save plot data --
    flux_barplot_save <- flux_barplot
    colnames(flux_barplot_save)[2:3] <- paste0("flux_mean_", c(unique(metadata_sub_group)[1],unique(metadata_sub_group)[2]))
    write.table(flux_barplot_save, paste0("03.scFEA_visualization/flux_", name, cluster, ".tsv"), sep="\t", row.names=F, quote=F)
    # set label for outliers --
    flux_barplot$label = trunc(flux_barplot$flux_change*100)/100   # two digits
    flux_barplot$flux_change[flux_barplot$flux_change>=2] = 2
    # plot --
    flux_barplot %>% filter(flux_change >0) %>%
        ggplot(aes(x=module_name, y=flux_change, fill=supermodule_name))+
            geom_bar(stat="identity", position="dodge", width = 0.7, color = "black")+
            geom_hline(yintercept = 1, color="red", linetype = 2)+
            geom_text(aes(y=abs(flux_change) + 0.03, label = label), size=1.5)+
            ylab(paste("flux_change:", unique(metadata_sub_group)[1], "vs", unique(metadata_sub_group)[2]))+
            coord_cartesian(ylim = c(0, 2.3)) +
            theme_classic()+theme(legend.position = "none", axis.title.x = element_blank(),
                                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 6),
                                  plot.title=element_text(size = 9))
        ggsave(paste0("03.scFEA_visualization/flux_barplot_", name, cluster, ".pdf"),device = "pdf", width = 20, height = 10)

    ## balance --
    balance_barplot <- data.frame(metabolite_name = factor(row.names(balance),levels = row.names(balance)),
                                  balance_g1 = apply(balance[,g1], 1, mean),
                                  balance_g2 = apply(balance[,g2], 1, mean),
                                  balance_change = apply(balance[,g1], 1, mean) / apply(balance[,g2], 1, mean))
    # set sign of balance  --> fill
    balance_barplot$direction <- rep("all negative",nrow(balance_barplot)) 
    balance_barplot$direction[balance_barplot$balance_g1 >0 & balance_barplot$ balance_g2 >0] = "all positive"
    balance_barplot$direction[balance_barplot$balance_g1 >0 & balance_barplot$ balance_g2 <0] = paste(unique(metadata_sub_group)[1],"positive")
    balance_barplot$direction[balance_barplot$balance_g1 <0 & balance_barplot$ balance_g2 >0] = paste(unique(metadata_sub_group)[2],"positive")                                
    balance_barplot$direction <- factor(balance_barplot$direction,levels=c("all positive", "all negative", paste(unique(metadata_sub_group)[1],"positive") ,paste(unique(metadata_sub_group)[2],"positive")))    
    # save plot data --
    balance_barplot_save <- balance_barplot
    colnames(balance_barplot_save)[2:3] <- paste0("balance_mean_", c(unique(metadata_sub_group)[1], unique(metadata_sub_group)[2]))
    write.table(balance_barplot_save,paste0("03.scFEA_visualization/balance_", name, cluster, ".tsv"), sep="\t", row.names=F, quote=F)
    # set labels --    
    balance_barplot$label = trunc(balance_barplot$balance_change*100)/100           
    balance_barplot$balance_change[abs(balance_barplot$balance_change)>=2] = 2  # set max 2      
    # plot --
    balance_barplot %>%
        ggplot(aes(x=metabolite_name, y=abs(balance_change), fill=direction))+
            geom_bar(stat="identity",position="dodge", width=0.7, alpha=0.7, color="black")+
            scale_fill_manual(values=c("#FF3030","#1E90FF","#FFA07A","#B0E2FF"))+
            geom_hline(yintercept=1, color="red", linetype=2)+
            geom_text(aes(y=abs(balance_change) + 0.03, label=label), size=1.5, color="black")+
            ylab(paste("balance_change:",unique(metadata_sub_group)[1],"vs",unique(metadata_sub_group)[2]))+
            coord_cartesian(ylim=c(0, 2.3)) +
            theme_classic()+theme(axis.title.x = element_blank(),
                                  axis.text.x = element_text(angle=90, hjust=1, vjust=0, size=6),
                                  plot.title=element_text(size=9))
    ggsave(paste0("03.scFEA_visualization/balance_barplot_", name, cluster,".pdf"),device="pdf", width=12, height=5)
}

function_scFEA_visualization <- function(cluster){
    cluster <- gsub(" ", "-", cluster)
    flux <- t(read.csv(paste0("02.scFEA_predication/", cluster, "_module_flux.csv"), row.names=1, check.names=F))  # results of scFEA
    balance <- t(read.csv(paste0("02.scFEA_predication/", cluster, "_balance.csv"), row.names=1, check.names=F))
    metadata_sub <- rds@meta.data[colnames(flux),]
    metadata_sub <- metadata_sub[order(metadata_sub[[group]]),]  # order by group for heatmap (deprecated)    
    flux <-  flux[,row.names(metadata_sub)] 
    row.names(flux) <- module_info[row.names(flux), 'Name']
    flux <- flux[c(1:4, 6, 5, 7:168),]  # swap M5 and M6 
    flux <- flux[apply(apply(flux, 1, duplicated), 2, sum)+1 !=nrow(metadata_sub),]  # rm module which all value equal
    balance <- balance[intersect(c("Glucose", word(module_info$Name, 2, sep=fixed("-> "))), row.names(balance)),]
    balance <- balance[c(1:5, 7, 6, 8:70),]  # swap Lactate and Acetyl-CoA
    # barplot --
    barplot1 <- function_barplot(model="split", flux, balance, metadata_sub, cluster)  # all data flux and balance
    cat(paste0(cluster, " done.\n"))
}
if(model == "split")    # barplot for scFEA using each cell type data
    visualization1 <- lapply(unique(rds@meta.data[[cell_type]]), function_scFEA_visualization)
if(model == "all")      # barplot for all data, 
    visualization1 <- function_scFEA_visualization("all")   

# subset all data fot visualization --
function_scFEA_visualization_all_sub <- function(cluster){
    cluster <- gsub(" ", "-", cluster)
    flux <- t(read.csv("02.scFEA_predication/all_module_flux.csv", row.names=1, check.names=F)) 
    balance <- t(read.csv("02.scFEA_predication/all_balance.csv", row.names=1, check.names=F))
    metadata_sub <- rds@meta.data[colnames(flux),]
    metadata_sub <- metadata_sub[order(metadata_sub[[group]]),]  

    flux <-  flux[,row.names(metadata_sub)]
    row.names(flux) <- module_info[row.names(flux),'Name']
    flux <- flux[c(1:4, 6, 5, 7:168),]  # swap M5 and M6 
    flux <- flux[apply(apply(flux, 1, duplicated), 2, sum)+1 !=nrow(metadata_sub),]  # rm module which all value equal
    balance <- balance[intersect(c("Glucose",word(module_info$Name, 2, sep=fixed("-> "))), row.names(balance)),]
    balance <- balance[c(1:5, 7, 6, 8:70),]  # swap Lactate and Acetyl-CoA
    # barplot --
    barplot2 <- function_barplot(model, flux, balance, metadata_sub, cluster)
    cat(paste0(cluster, " done.\n"))
}
if(model=="all"){
    cat("Split flux and balance by cell type, start visuslization...\n")
    visualization2 <- lapply(unique(rds@meta.data[[cell_type]]), function_scFEA_visualization_all_sub)  # all each cluster
}

# print run time --
runningtime <- proc.time()- timestart
cat("----\n")
cat(paste0("Running time: ", as.numeric(runningtime[1])," secs. \n"))  