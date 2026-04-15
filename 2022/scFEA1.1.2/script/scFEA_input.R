## code to generate input of scFEA ----

cat("scFEA start running...\n")
#timestart <- proc.time()

# library packages --
suppressMessages(library(Seurat))
suppressMessages(library(argparser))
suppressMessages(library(tidyverse))

# args ----
args <- arg_parser('')
args <- add_argument(args,"--rds", help="seurat rds file")
args <- add_argument(args,"--exp", help="exp matrix file")
args <- add_argument(args,"--model", help="model,split or all")
args <- add_argument(args,"--cell_type", help="name of variable about cell type")
args <- add_argument(args,"--group", help="name of variable about group")
args <- add_argument(args,"--species", help="species")
args <- add_argument(args,"--n_threads", help="num of threads to run scFEA")
args <- add_argument(args,"--sc_imputation", help="sc_imputation")
args <- add_argument(args,"--moduleGene_file", help="moduleGene_file")
args <- add_argument(args,"--stoichiometry_matrix", help="stoichiometry_matrix")
args <- add_argument(args,"--cName_file", help="cName_file")
args <- add_argument(args,"--outdir", help="outdir")
args <- parse_args(args)


MAX_samples = 10^4
set.seed(16)   # set seed for downsample

# read rds --
if(args$rds != "None"){
    rds = readRDS(args$rds)
}
if(args$exp != "None"){
    exp = read.table(args$exp,row.names=1,header=T,sep="\t")
}
model <- args$model
cell_type <- args$cell_type
group <- args$group
species <- args$species
n_threads <- args$n_threads
sc_imputation <- args$sc_imputation
moduleGene_file <- args$moduleGene_file
stoichiometry_matrix <- args$stoichiometry_matrix
cName_file <- args$cName_file
outdir <- args$outdir

# mkdir --
if(!dir.exists(outdir))
    dir.create(outdir)
if(!dir.exists(paste0(outdir,"/01.scFEA_input")))
    dir.create(paste0(outdir,"/01.scFEA_input"))
if(!dir.exists(paste0(outdir,"/02.scFEA_predication/")))
    dir.create(paste0(outdir,"/02.scFEA_predication/"))
if(!dir.exists(paste0(outdir,"/03.scFEA_visualization/")))
    dir.create(paste0(outdir,"/03.scFEA_visualization/"))

# set 3 files --
if(moduleGene_file == "None" | stoichiometry_matrix == "None" | cName_file == "None"){
    if(species == "human"){
        moduleGene_file="module_gene_m168.csv"
        stoichiometry_matrix="cmMat_c70_m168.csv"
        cName_file="cName_c70_m168.csv"  
    }
    if(species == "mouse"){
        moduleGene_file="module_gene_complete_mouse_m168.csv"
        stoichiometry_matrix="cmMat_complete_mouse_c70_m168.csv"
        cName_file="cName_complete_mouse_c70_m168.csv" 
    }
}

#
if(args$rds != "None"){
    print("rds")
    cell_type <- rds@meta.data[,cell_type]
    n_cell_types = length(unique(cell_type))
    cell_type_large <- names(table(cell_type)[(table(cell_type) > MAX_samples)])

    # downsample --
    function_large_downsample <- function(cluster){
        metadata_sub <- rds@meta.data[cell_type == cluster,]  
        group_downsample <- table(metadata_sub[,group])/nrow(metadata_sub)* MAX_samples
        downsample_cells <- c()
        for( i in names(group_downsample)){
            group_cells <- row.names(metadata_sub)[metadata_sub[,group]==i]
            downsample_cells <- c(downsample_cells, sample(group_cells)[1:as.numeric(floor(group_downsample[i]))])
        }  
        return(downsample_cells)
    }
    # save data --
    function_save_by_cell_type <- function(cluster){
        if(cluster %in% cell_type_large){
            cat(paste0(cluster,"> ", MAX_samples, " cells, start downsampling\n"))
            data <- rds@assays$RNA@data[,match(function_large_downsample(cluster), colnames(rds@assays$RNA@data))]
        }
        else
            data <- rds@assays$RNA@data[,cell_type == cluster]
        write.csv(data, paste0(outdir,"/01.scFEA_input/input_", gsub(" ", "-", cluster), ".csv"), quote=F)
        cat(paste0(cluster, " done!\n"))
    }
    if(model == "split"){
        save_data <- lapply(unique(cell_type), function_save_by_cell_type)
        mapfile <- data.frame(test_file=paste0(outdir,"/01.scFEA_input/input_", gsub(" ", "-", unique(cell_type)), ".csv"),
                              species = rep(species, n_cell_types),
                              result_prefix = gsub(" ", "-", unique(cell_type)),
                              moduleGene_file = rep(moduleGene_file, n_cell_types),
                              stoichiometry_matrix = rep(stoichiometry_matrix, n_cell_types),
                              cName_file = rep(cName_file, n_cell_types),
                              n_threads = rep(n_threads, n_cell_types),
                              sc_imputation = rep(sc_imputation, n_cell_types))
        write.csv(mapfile, file=paste0(outdir, "/01.scFEA_input/mapfile.csv"), quote=F, row.names=F)
    }
    if(model == "all"){
        write.csv(rds@assays$RNA@data, paste0(outdir, "/01.scFEA_input/input_all.csv"), quote=F)
        mapfile <- data.frame(test_file=paste0(outdir, "/01.scFEA_input/input_all.csv"),
                              species = species,
                              result_prefix = "all",
                              moduleGene_file = moduleGene_file,
                              stoichiometry_matrix = stoichiometry_matrix,
                              cName_file = cName_file,
                              n_threads = n_threads,
                              sc_imputation = sc_imputation)
        write.csv(mapfile, file=paste0(outdir, "/01.scFEA_input/mapfile.csv"), quote=F, row.names=F)
    }
}

if(args$exp != "None"){
    print("exp matrix")
    write.csv(exp, paste0(outdir, "/01.scFEA_input/input_all.csv"), quote=F)
    mapfile <- data.frame(test_file=paste0(outdir, "/01.scFEA_input/input_all.csv"),
                          species = species,
                          result_prefix = "all",
                          moduleGene_file = moduleGene_file,
                          stoichiometry_matrix = stoichiometry_matrix,
                          cName_file = cName_file,
                          n_threads = n_threads,
                          sc_imputation = sc_imputation)
    write.csv(mapfile, file=paste0(outdir, "/01.scFEA_input/mapfile.csv"), quote=F, row.names=F)
}

# print run time --
#runningtime <- proc.time()- timestart
#cat("----\n")
#cat(paste0("Running time: ", as.numeric(runningtime[1])," secs. \n"))  