suppressMessages({
    library(Seurat)
    library(argparser)
})

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--virus_mtx", help = "virus mtx, split by ,")
argv <- add_argument(argv, "--prefix", help = "zl sample name, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir. Deafalut: ./")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

if(is.na(argv$outdir)){
    argv$outdir ="./"
}
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir)
}

virus_mtx <- unlist(strsplit(argv$virus_mtx, split = ","))
prefix <- unlist(strsplit(argv$prefix, split = ","))

print(paste0("virus mtx: ", virus_mtx))
print(paste0("prefix: ", prefix))
print("--------")

# read rds -- 
data_seurat <- readRDS(argv$rds)

function_read_virus_mtx <- function(i){
    data <- Read10X(virus_mtx[i])
    colnames(data) <- paste0(prefix[i], "_", colnames(data))
    return(data)
}
mtx <- do.call(cbind, lapply(1:length(virus_mtx), function_read_virus_mtx))

if(!identical(colnames(mtx),colnames(data_seurat))){
    print("The cell id between rds and virus mtx is different. Try subset ...")
    
    if(length(intersect(colnames(mtx), colnames(data_seurat))) == ncol(data_seurat)){
        print("The cell id between rds and (subset) virus mtx is same.")
        mtx <- mtx[, colnames(data_seurat)]  # match
        for (i in row.names(mtx)){
            data_seurat@meta.data[i] <- mtx[i,]
        }
        saveRDS(data_seurat, paste0(argv$outdir, "/", argv$name ,".rds"))
        print("Done.")
    }else{
        print("The cell id between rds and virus mtx is still different. Please check input.")
    }
}else{
    print("The cell id between rds and virus mtx is same.")
    for (i in row.names(mtx)){
        data_seurat@meta.data[i] <- mtx[i,]
    }
    saveRDS(data_seurat, paste0(argv$outdir, "/", argv$name ,".rds"))
    print("Done.")
}