options(warn = -1)    # off warnings 
suppressMessages({
    library(argparser)
    library(tidyverse)})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--network", help = "network file")
argv <- add_argument(argv, "--node", help = "node file")
argv <- add_argument(argv, "--regulons", help = "regulons")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

#
network <- read.table(argv$network, sep="\t", header = T) 
network_subset <- network[ network$type %in% paste0('interacts with ', unlist(strsplit(argv$regulon, split = ","))),]
write.table(network_subset, paste0(argv$outdir, "/cytoscape_regulon_subset_network.txt"), sep="\t", quote = F, row.names = F)

node <- network_subset[ !duplicated(paste0(network_subset$target, network_subset$type)),]
multi_target <- unique(node$target[ duplicated(node$target)])   # multi tf target same gene
node$gene_type = 'target'
colnames(node)[2]='gene'

# tf behind, 
node_df <- rbind(node[,c('gene','gene_type','type')],
                 data.frame(gene = unique(node$tf),
                            gene_type = 'tf',
                            type = paste0("interacts with ",unique(node$tf))))
node_df$type[ node_df$gene %in% multi_target] = 'multi'
write.table(node_df, paste0(argv$outdir, "/cytoscape_regulon_subset_node.txt"),sep="\t",quote = F,row.names = F)

print('regulons to cytoscape done.')