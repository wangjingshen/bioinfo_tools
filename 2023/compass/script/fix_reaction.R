options(warn = -1)
suppressMessages({
    library(argparser)
    library(tidyverse)
    library(data.table)
})
options(warn=1)

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--reactions", help = "compass reactions")
argv <- add_argument(argv, "--metadata", help = "compass reaction metadata")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

argv$metadata <- ifelse(is.na(argv$metadata), "/SGRNJ06/randd/USER/wangjingshen/script/compass/data/reaction_metadata.csv", argv$metadata)
argv$outdir <- ifelse(is.na(argv$outdir), "./", argv$outdir)
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir)
}
argv$name <- ifelse(is.na(argv$name), "", paste0("_", argv$name))

# read --
res = read.table(argv$reactions, sep = "\t", header = T, row.names = 1)
ref <- fread(argv$metadata, encoding = 'UTF-8')

# function from compass
function_get_reaction_consistencies <- function(compass_reaction_penalties, min_range=1e-3){
    df = -log(compass_reaction_penalties + 1)
    df = df[apply(res,1,max) - apply(res,1,min) >= min_range,]
    df = df - min(df)
    return(df)
}
res <- function_get_reaction_consistencies(res)

res_ref = ref[match(sapply(row.names(res), function(x){substring(x, 1, nchar(x) - 4)}), ref$reaction_no_direction),]
row.names(res) <- paste0(row.names(res), "___", res_ref$subsystem, "___", res_ref$reaction_name)

# out
write.table(res, paste0(argv$outdir, "/reaction_consistencies", argv$name, ".tsv"), quote = F, sep="\t")