suppressWarnings(suppressMessages({
    library(argparser)
    library(ArchR)
    library(SingleR)
}))
set.seed(1)
addArchRThreads(threads = 4)


argv <- arg_parser('')
argv <- add_argument(argv,"--fragments", help="fragments file")
argv <- add_argument(argv,"--sname", help="sname")
argv <- add_argument(argv,"--rds", help="rds")
argv <- add_argument(argv,"--species", help="species, Default: human")
argv <- add_argument(argv,"--minTSS", help="minTSS, Default: 4")
argv <- add_argument(argv,"--minFrags", help="minFrags, Default: 1000")
argv <- add_argument(argv,"--outdir", help="outdir, Default: outdir")
argv <- parse_args(argv)

species <- ifelse(is.na(argv$species), "human", argv$species)
minTSS <- ifelse(is.na(argv$minTSS), 4, as.numeric(argv$minTSS))
minFrags <- ifelse(is.na(argv$minFrags), 1000, as.numeric(argv$minFrags))
outdir <- ifelse(is.na(argv$outdir), "outdir", as.numeric(argv$outdir))
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

if(is.na(argv$rds)){
    fragments <- unlist(strsplit(argv$fragments, split = ","))
    sname <- unlist(strsplit(argv$sname, split = ","))
}else{
    rds <- argv$rds
}

if(species == "human"){
    addArchRGenome("hg38")
}else{
    addArchRGenome("mm10")
}

function_singleR <- function(data, species){
    if(species == "human"){
        ref = readRDS("/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/HumanPrimaryCellAtlasData.rds")
    }
    if(species == "mouse"){
        ref = readRDS("/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/MouseRNAseqData.rds")
    }
    data_gene_score <- getMatrixFromProject(ArchRProj = data, useMatrix = "GeneScoreMatrix", useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
    # identical(colnames(data_gene_score),row.names(data@cellColData))
    gene_score <- as.matrix(data_gene_score@assays@data$GeneScoreMatrix)
    row.names(gene_score) <- rowData(data_gene_score)$name
    anno <- SingleR(test = gene_score, ref = ref, 
                    clusters = unlist(data$Clusters ),
                    assay.type.test=1, labels = ref$label.main)
    cell_type_singleR <- plyr::mapvalues(x = data$Clusters, from = row.names(anno), to = anno$labels)
    return(cell_type_singleR)
}

if(is.na(argv$rds)){
    ArrowFiles <- createArrowFiles(inputFiles = fragments, sampleNames = sname, minTSS = 4, minFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)
    doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)
    proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "result", copyArrows = TRUE)
    # proj <- filterDoublets(ArchRProj = proj)
    proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
    proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
    proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
    # getAvailableMatrices(proj)

    p1 <- plotEmbedding(ArchRProj = data, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
    p2 <- plotEmbedding(ArchRProj = data, colorBy = "cellColData", name = "cell_type_singleR", embedding = "UMAP")
    ggAlignPlots(p1, p2, type = "h")
    plotPDF(p1, p2, name = paste0(outdir, "/anno.pdf"), ArchRProj = data, addDOC = FALSE, width = 8, height = 5)
}else{
    data <- readRDS(argv$rds)
    data@projectMetadata$outputDirectory = "outdir"
    data$cell_type_singleR <- function_singleR(data, "human")

    p1 <- plotEmbedding(ArchRProj = data, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
    p2 <- plotEmbedding(ArchRProj = data, colorBy = "cellColData", name = "cell_type_singleR", embedding = "UMAP")
    ggAlignPlots(p1, p2, type = "h")
    plotPDF(p1, p2, name = "anno.pdf", ArchRProj = data, addDOC = FALSE, width = 8, height = 5)
}
