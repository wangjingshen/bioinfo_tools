options(warn = -1) 
suppressMessages({
    library(tidyverse)
    library(argparser)
})
options(warn = 1) 


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--mapfile", help = "mapfile, split by , ")
argv <- add_argument(argv, "--subset", help = "subset sample, split by , ")
argv <- add_argument(argv, "--mode", help = "mode, zl or fj")
argv <- add_argument(argv, "--time", help = "time")
argv <- parse_args(argv)

# zl
mapfile <- unlist(strsplit(argv$mapfile, split = ","))
mapfile_df <- do.call(rbind, lapply(mapfile,function(x){   read.table(x,sep="\t")   }))
if(argv$subset != "None"){ 
    mapfile_df <- mapfile_df[mapfile_df[,3] %in% unlist(strsplit(argv$subset, split = ",")), ]
}
mapfile_df[,3] <- gsub("HBVV[12]_", "", mapfile_df[,3])  # gsub HBVV1 HBVV2
name <- paste0(mapfile_df[,3], "---", word(mapfile_df[,2], 5, sep = fixed("/")))     # sample_id --- time
name[duplicated(name)] <- paste0(name[duplicated(name)], "_", 1:sum(duplicated(name)))

# sub time 
if(sum(grepl("RD20102301_DZH",name))!=0){
    cat(paste0("---\nExist non-time-name: ", name[grepl("RD20102301_DZH",name)],". Replace with ", argv$time, ".\n---\n"))
    name <- gsub("RD20102301_DZH", argv$time, name)
}


# new mapfile
if(argv$mode == "zl"){
    mapfile_df_copy = mapfile_df
    mapfile_df_copy[,2] <- paste0("fastq/ZL/", name)
    write.table(mapfile_df_copy, file = "data/zl_mapfile", sep="\t", row.names = F, col.names = F, quote = F)
}

if(argv$mode == "fj"){
    mapfile_df_copy = mapfile_df
    mapfile_df_copy[,4] <- gsub("HBVV[12]_", "", mapfile_df_copy[,4])  # gsub HBVV1 HBVV2
    mapfile_df_copy[,2] <- paste0("fastq/FJ/", name)
    mapfile_df_copy[,4] <- paste0(word(mapfile_df_copy[,4], -1, sep=fixed("/")), "_celescope_dir")    # zl celescope dir
    write.table(mapfile_df_copy, file = "data/fj_mapfile", sep="\t", row.names = F, col.names = F, quote = F)
}
