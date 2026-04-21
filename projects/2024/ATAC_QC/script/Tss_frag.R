suppressWarnings(suppressMessages({
    library (readr)
    library (argparser)
    library(ArchR)
    library(pheatmap)
    library (Seurat)
    library (tidyr)
    library (dplyr)
    library (stringr)
}))  


## 物种还可以是果蝇
argv <- arg_parser('')
argv <- add_argument (argv,"--fragment", help="the cell gene express file list, split by ,")
argv <- add_argument (argv,"--spname",help="the samples name list ,split by ,")
argv <- add_argument (argv,"--gname",help = "the group name list,split by ,")
argv <- add_argument (argv,"--info_table",help = "scATAC info",default = "F")
argv <- add_argument (argv,"--species",help="hg19 or hg38 or mm10 or mm9 or pig or Dmelanogaster",default = "hg38")
argv <- add_argument (argv,"--outdir",help="path of outdir")
argv <- add_argument (argv,"--filterTSS",help = "TSS filter threshold",default = 4)  # 4
argv <- add_argument (argv,"--filterFrags",help = "Fragment filter threshold",default = 1000)    # 1000

argv <- parse_args(argv)
info_table <- argv$info_table
outdir <- argv$outdir
dir.create (outdir,showWarnings=T)

set.seed (1111)
filterTSS <- argv$filterTSS
filterFrags <- argv$filterFrags
#obj <- as.numeric(argv$object)
# wd <- getwd()
## Setting a Genome and GeneAnnotation
if (argv$species == "hg38") {
    addArchRGenome("hg38")
} else if (argv$species == "hg19") {
    addArchRGenome("hg19")
} else if (argv$species == "mm10") {
    addArchRGenome("mm10")
} else if (argv$species == "mm9") {
    addArchRGenome("mm9")
}


# ### ignore the chromosome prefixes
# #addArchRChrPrefix(chrPrefix = FALSE)

# ## Setting default number of Parallel threads.
addArchRThreads(threads = 1)
print ("Creating Arrow Files")
#### Creating Arrow Files
### Per-cell Quality Control (TSS enrichment score greater than 4 and more than 1000 unique nuclear fragments)
if (info_table != "F") {
    info <- read.table (info_table,header = T,sep='\t')
    samplename <- info$spname
    inputFiles <- info$fragment
    groupname <- info$gname
    names (groupname) <- samplename
}else {
    inputFiles <- c()
    if (grepl(',',argv$fragment)) {
        inputFiles <- unlist(strsplit(argv$fragment,split=","))
    }else {
        inputFiles <- c(argv$fragment)
    }
    samplename <- c()
    if (grepl(",",argv$spname)) {
        samplename <- unlist(strsplit(argv$spname,split=','))
    }else {
        samplename <- c(argv$spname)
    }
    groupname <- c()
    if (grepl(",",argv$gname)) {
        groupname <- unlist(strsplit(argv$gname,split=','))
    }else {
        groupname <- c(argv$gname)
    }
    names (groupname) <- samplename
}


names (inputFiles) <- samplename
print (inputFiles)

ArrowFiles <- createArrowFiles(
inputFiles = inputFiles,
sampleNames = names(inputFiles),
filterTSS = filterTSS, #Dont set this too high because you can always increase later
filterFrags = filterFrags,
addTileMat = TRUE,
addGeneScoreMat = TRUE
)
ArrowFiles

## Inferring scATAC-seq Doublets (数据需不需要去doublet需要根据R2 确认)
### This is likely a good place to start but we encourage all users to inspect the pre- and post-doublet removal data to understand how doublet removal is affecting the cells
### ArchR reports the R2 value for the UMAP projection for each Arrow file.
### If these R2 values are much lower (i.e. less than 0.9), skipping doublet prediction.
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

### Creating an ArchRProject
projHeme1 <- ArchRProject(
ArrowFiles = ArrowFiles,
outputDirectory = outdir,
copyArrows = TRUE  #This is recommened so that if you modify the Arrow files you have an original copy for later usage.##看一下修改为F结果是怎么样的
)
projHeme1
#### which data matrices are available within the ArchRProject
getAvailableMatrices(projHeme1)
###########画图
### Plotting QC metrics
#### log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)
dir.create (file.path(outdir,"Plots"),showWarnings=T)
dir.create (file.path(outdir,"Plots","0.QCmetrics"),showWarnings=T)
dip <- file.path (outdir,"Plots","0.QCmetrics")
plotPDF(p, name = file.path("0.QCmetrics","TSS-vs-Frags.pdf"), ArchRProj = projHeme1, addDOC = FALSE)
png(file.path(dip,"TSS-vs-Frags.png"))
p
dev.off()


# plot TSS
#p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
#p2
#dip <- file.path (outdir,"Plots","0.QCmetrics")
#plotPDF(p2, name = file.path("0.QCmetrics","TSS.pdf"), ArchRProj = projHeme1, addDOC = FALSE)
#png(file.path(dip,"TSS.png"))
#p2
#dev.off()