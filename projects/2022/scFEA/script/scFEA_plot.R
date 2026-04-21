## script for a simple example of downstream analysis of scFEA ---

options(stringsAsFactors = FALSE)

# library packages --
suppressMessages(library(tidyverse))
suppressMessages(library(pheatmap))
suppressMessages(library(argparser))

# args --
args <- arg_parser('')
args <- add_argument(args,"--flux_module", help="modlue flux from scFEA")
args <- add_argument(args,"--balance_metabolite", help="metabolite balance from scFEA")
args <- add_argument(args,"--metadata", help="metadata including group info")
args <- add_argument(args,"--superModule_sub", default = "../data/superModule_sub.csv", help="interest supermodules")
args <- add_argument(args,"--metabolite_sub", default = "../data/metabolite_sub.csv", help="interest metabolie")
args <- add_argument(args,"--output_path", default = "./", help="path to save results")
args <- parse_args(args)


## read data ----
flux_module <- t(read.csv(args$flux_module,row.names = 1, check.names = F))
balance_metabolite <- t(read.csv(args$balance_metabolite,row.names = 1, check.names = F))

ref_module <- read.csv("/SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA/data/module_info_add_supermodule.csv",row.names = 1)
metadata <- read.csv(args$metadata,row.names = 1)
g1 <- metadata$cell[metadata$group == unique(metadata$group)[1]]
g2 <- metadata$cell[metadata$group == unique(metadata$group)[2]]


## plot ----
# heatmap ---
superModule_sub <- read.csv(args$superModule_sub,header=F)
heatmap_data <- t(flux_module[row.names(flux_module) %in% row.names(ref_module)[ref_module$Supermodule_name %in% superModule_sub],])
colnames(heatmap_data) <- ref_module$Module_name[ref_module$Supermodule_name %in% superModule_sub]
annotation_row <- data.frame(row.names = metadata$cell,
                             group = metadata$group)
pdf(paste0("scFEA_plot/","heatmap_fulx_between_two_groups.pdf"))
pheatmap(heatmap_data,show_rownames = F,cluster_rows = F,cluster_cols = F,scale = "row",angle_col = 90,annotation_row = annotation_row)
dev.off()

# flux change barplot --
data.frame(module_name = ref_module[row.names(flux_module),]$Module_name,
           supermodule_name = ref_module[row.names(flux_module),]$Supermodule_name,
           flux_change = apply(flux_module[,g1], 1, sum) / apply(flux_module[,g2], 1, sum)) %>%
  filter(supermodule_name %in% superModule_sub) %>%
  ggplot(aes(x=module_name,y=flux_change))+
  geom_bar(stat="identity",position="dodge",fill="#00BFFF",width = 0.7,color = "black")+
  geom_hline(yintercept = 1, color="red")+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
ggsave(paste0(args$output_path,"barplot_fulx_change_between_two_groups.pdf"),device = "pdf")

# metabolic change barplot--
metabolite_sub <- read.csv(args$metabolite_sub, header=F)
data.frame(metabolite_name = row.names(balance_metabolite) ,
           flux_change = apply(balance_metabolite[,g1], 1, sum) / apply(balance_metabolite[,g2], 1, sum)) %>%
  filter(metabolite_name %in% metabolite_sub) %>%
  ggplot(aes(x=metabolite_name,y=flux_change))+
  geom_bar(stat="identity",position="dodge",fill="#00BFFF",width = 0.7,color = "black")+
  geom_hline(yintercept = 1, color="red")+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1),
                        axis.title.x = element_blank())
ggsave(paste0(args$output_path,"barplot_metabolic_change_between_two_groups.pdf"),device = "pdf")
