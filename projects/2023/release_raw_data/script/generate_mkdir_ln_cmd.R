options(warn = -1) 
suppressMessages({
    library(tidyverse)
    library(argparser)
})
options(warn = 1) 


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--mapfile", help = "mapfile, split by , ")
argv <- add_argument(argv, "--project_id", help = "project id")
argv <- add_argument(argv, "--subset", help = "subset sample, split by , ")   # for subset release
argv <- add_argument(argv, "--mode", help = "mode, zl or fj")
argv <- add_argument(argv, "--rm_version", help = "T or F")
argv <- add_argument(argv, "--outdir", help = "outdir")
argv <- add_argument(argv, "--generate_mapfile", help = "generate mapfile. Default: no")
argv <- add_argument(argv, "--time", help = "time")
argv <- parse_args(argv)

argv$generate_mapfile <- ifelse(is.na(argv$generate_mapfile), 'no', argv$generate_mapfile)

# cmd dir
if(!dir.exists("cmd/")){
    dir.create("cmd/")
}
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = T)
}
cat(paste0("mode: ", argv$mode,"\n"))

mapfile <- unlist(strsplit(argv$mapfile, split = ","))
cat(c("mapfile:\n", paste0(mapfile,"\n")))
#cat(paste0("mapfile: ", mapfile,"\n"))
if(argv$subset != "None"){
    cat(paste0("subset: ", unlist(strsplit(argv$subset, split = ",")),"\n"))
}

mapfile_df <- do.call(rbind, lapply(mapfile,function(x){   read.table(x,sep="")   }))
print(mapfile_df)
if(argv$subset != "None"){
    mapfile_df <- mapfile_df[mapfile_df[,3] %in% unlist(strsplit(argv$subset, split = ",")), ]
}
if(argv$rm_version == "T"){
    mapfile_df[,3] <- gsub("HBVV[12]_", "", mapfile_df[,3])  # gsub HBVV1 HBVV2
}
name <- paste0(mapfile_df[,3], "---", word(word(mapfile_df[,2], 1, sep = fixed(paste0("/", argv$project_id))), -1, sep=fixed("/")))     # sample_id-time
name[duplicated(name)] <- paste0(name[duplicated(name)], "_", 1:sum(duplicated(name)))
print(name)

argv$rm_version <- ifelse(is.na(argv$rm_version), "F", argv$rm_version)
print(argv$rm_version)


# sub time
if(!is.na(argv$time)){
    if(sum(grepl("RD20102301_DZH",name))!=0){
        cat(paste0("---\nExist non-time-name: ", name[grepl("RD20102301_DZH",name)],". Replace with ", argv$time, ".\n---\n"))
        name <- gsub("RD20102301_DZH", argv$time, name)
    }
}

if(argv$mode == "zl"){
    # cmd
    mkdir_cmd <- paste0("mkdir -p ", argv$outdir, "/fastq/ZL/", name)
    ln_cmd <- paste0("ln -s ", mapfile_df[,2], "/", mapfile_df[,1], "*fastq.gz ", argv$outdir, "/fastq/ZL/", name)
    write.table(mkdir_cmd, "cmd/zl_mkdir.sh", sep="\n", row.names = F, col.names = F, quote = F)
    write.table(ln_cmd, "cmd/zl_ln.sh", sep="\n", row.names = F, col.names = F, quote = F)
    # mapfile
    if(argv$generate_mapfile != 'no'){
        mapfile_df_copy = mapfile_df
        mapfile_df_copy[,2] <- paste0(argv$outdir, "fastq/ZL/", name)
        write.table(mapfile_df_copy, file = paste0(argv$outdir, "/zl_mapfile"), sep="\t", row.names = F, col.names = F, quote = F)
        zip_mapfile_cmd <- paste0('zip -m ', argv$outdir, "/zl_mapfile.zip ", argv$outdir, "/zl_mapfile")
        write.table(zip_mapfile_cmd, "cmd/zl_mapfile_zip.sh", sep="\n", row.names = F, col.names = F, quote = F)
    }
}
if(argv$mode == "fj"){
    # cmd
    mkdir_cmd <- paste0("mkdir -p ", argv$outdir, "/fastq/FJ/", name)
    ln_cmd <- paste0("ln -s ", mapfile_df[,2], "/", mapfile_df[,1], "*fastq.gz ", argv$outdir, "/fastq/FJ/", name)
    write.table(mkdir_cmd, "cmd/fj_mkdir.sh", sep="\n", row.names = F, col.names = F, quote = F)
    write.table(ln_cmd, "cmd/fj_ln.sh", sep="\n", row.names = F, col.names = F, quote = F)
    # mapfile
    if(argv$generate_mapfile != 'no'){
        mapfile_df_copy = mapfile_df
        mapfile_df_copy[,2] <- paste0(paste0(argv$outdir, "/fastq/FJ/"), name)
        mapfile_df_copy[,4] <- paste0(word(mapfile_df_copy[,4], -1, sep=fixed("/")), "_celescope_dir")    # zl celescope dir
        write.table(mapfile_df_copy, file = paste0(argv$outdir, "/fj_mapfile"), sep="\t", row.names = F, col.names = F, quote = F)
        zip_mapfile_cmd <- paste0('zip -m ', argv$outdir, "/fj_mapfile.zip ", argv$outdir, "/fj_mapfile")
        write.table(zip_mapfile_cmd, "cmd/fj_mapfile_zip.sh", sep="\n", row.names = F, col.names = F, quote = F)
    }
}