
options(warn = -1)    # off warnings 
suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
})
options(warn = 1)     # 

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--mix", help = "mix")
argv <- add_argument(argv, "--zl_matrix", help = "zl matrix")
argv <- add_argument(argv, "--outdir", help = "outdir, Default: outdir/")
argv <- parse_args(argv)

if(!is.na(argv$outdir)){
    argv$outdir = "outdir/"
}
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir)
}

mix <- read.table(argv$mix,sep="\t",header = T)
count = Read10X(argv$zl_matrix)
mix[1,1] <- colnames(count)[1]   # fix bug of doublet mix
cat(paste0("mix cell_id and matrix barcode: ",identical(mix[,1], colnames(count)),"\n"))
print(table(mix$Species))
hs <- count[, colnames(count) %in% mix$cell_id[ mix$Species=="human"]]
mm <- count[, colnames(count) %in% mix$cell_id[ mix$Species=="mouse"]]
cat(paste0("hs matrix ncol: ", ncol(hs),"\n"))
cat(paste0("mm matrix ncol: ", ncol(mm),"\n"))
res_df = data.frame(variable = c("double_percent", "mean_UMI_hs","mean_genes_hs","mean_UMI_mm","mean_genes_mm"),
                    values = c(paste0(round(sum(mix$Species == "multiplet")/ncol(count)*100,2)," %"),
                               round(mean(colSums(hs))),
                               round(mean(colSums(hs>0))),
                               round(mean(colSums(mm))),
                               round(mean(colSums(mm>0)))))
write.table(res_df, paste0(argv$outdir, "/doublet_stat.txt"), sep="\t",quote=F, row.names=F ,col.names=F)

cat("----------------\ndoublet stat done.")