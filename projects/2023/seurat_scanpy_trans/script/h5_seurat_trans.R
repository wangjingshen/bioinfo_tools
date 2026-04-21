suppressMessages({
    library(Seurat)
    library(dior)
    library(argparser)
})

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--h5", help = "h5")
argv <- add_argument(argv, "--input_type", help = "input_type, rds or h5")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

# seurat to h5 (to scanpy) --
if(argv$input_type == 'rds'){
    data_seurat <- readRDS(argv$rds)
    write_h5(data_seurat, file = paste0(argv$outdir, '/', argv$name, '.h5'),
        object.type='seurat',assay.name='RNA',save.graphs=TRUE, save.scale=FALSE)
}

# (scanpy to) h5 to seurat --
if(argv$input_type == 'h5'){
    data_seurat <- read_h5(file = argv$h5,assay.name = 'RNA', target.object = 'seurat')
    saveRDS(data_seurat, paste0(argv$outdir, '/', argv$name, ".rds"))
}
