suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

argv <- arg_parser('')
argv <- add_argument(argv, "--mtx", help = "path of celescope space mtx")
argv <- add_argument(argv, "--pos", help = "tissue_positions")
argv <- add_argument(argv, "--spname", help = "spname")
argv <- add_argument(argv, "--swap_pos", help = "swap_pos")
argv <- parse_args(argv)

mtx <- argv$mtx
pos <- argv$pos
spname <- argv$spname
swap_pos <- argv$swap_pos

# make cnts
data <- Read10X(mtx)
write.table(t(as.matrix(data)), str_glue("{spname}/cnts.tsv"), sep="\t", quote = F)

# make pos
pos <- read.table(pos, sep=",", header = 0)
if(swap_pos %in% c("T","True", "TRUE")){
    pos <- pos[,c(1, 6, 5)]
}else{
    pos <- pos[,c(1, 5, 6)]
}
colnames(pos) <- c("spot","x","y")
pos <- pos[ match(colnames(data), pos$spot),]
write.table(pos, str_glue("{spname}/locs-raw.tsv"), sep="\t", quote = F, row.names = F)

print("generate cnts locs done.")