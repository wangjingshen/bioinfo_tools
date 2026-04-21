options(warn = -1)    # off warnings 
suppressMessages({
    library(argparser)
    library(tidyverse)})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--reg", help = "reg csv from pyscenic")
argv <- add_argument(argv, "--subset", help = "subset sample, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

# 
reg <- read.table(argv$reg, sep=',')
colnames(reg) <- cbind(reg[3,1:2], reg[2,3:10])
reg <- reg[-(1:3),]
function_get_regulon <- function(i){
    target_list = gsub(" ","",unlist(str_split(gsub('\'|\\[|\\]|\\(|\\)', '', reg$TargetGenes[i]), ',')))
    target_df = data.frame(target = target_list[seq(1,length(target_list),2)],
                           score = target_list[seq(2,length(target_list),2)])
    target_df$tf = reg$TF[i]
    target_df$type <- paste0('interacts with ',target_df$tf)
    return(target_df[,c(3,1,2,4)])
}
regulon_df <- do.call(rbind, lapply(1:nrow(reg),function_get_regulon))
regulon_df <- regulon_df[!duplicated(paste0(regulon_df$tf,regulon_df$target)),]
# network + edge(type)
if(!is.na(argv$subset)){
    regulon_df <- regulon_df[ regulon_df$tf %in% unlist(strsplit(argv$subset, split = ",")),]
}
write.table(regulon_df[,c(1,2,4)], paste0(argv$outdir,"/cytoscape_network.txt"), sep="\t", quote = F, row.names = F)

regulon_df_node <- regulon_df[ !duplicated(paste0(regulon_df$target, regulon_df$type)),c(1,2,4)]
multi_target <- unique(regulon_df_node$target[ duplicated(regulon_df_node$target)])   # multi tf target same gene
regulon_df_node$gene_type = 'target'
colnames(regulon_df_node)[2]='gene'

# tf behind, 
node_df <- rbind(regulon_df_node[,c(2,4,3)],
                 data.frame(gene = unique(regulon_df_node$tf),
                            gene_type = 'tf',
                            type = paste0("interacts with ",unique(regulon_df_node$tf))))
node_df$type[ node_df$gene %in% multi_target] = 'multi'
write.table(node_df, paste0(argv$outdir, "/cytoscape_node.txt"),sep="\t",quote = F,row.names = F)

print('regulons to cytoscape done.')