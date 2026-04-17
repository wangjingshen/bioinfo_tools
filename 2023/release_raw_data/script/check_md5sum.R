options(warn = -1)
suppressMessages({
    library(argparser)
    library(tidyverse)
})
options(warn = 1)

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--path", help = "path of md5sum file")
argv <- parse_args(argv)


if(file.exists("cmd/zl_ln.sh") & file.exists("cmd/fj_ln.sh")){
    ln_cmd <- rbind(read.table("cmd/zl_ln.sh", sep=" "),
                    read.table("cmd/fj_ln.sh", sep=" "))
}else if(file.exists("cmd/zl_ln.sh")){
    ln_cmd <- read.table("cmd/zl_ln.sh", sep=" ")
}else{
    ln_cmd <- read.table("cmd/fj_ln.sh", sep=" ")
}

function_read_raw_md5 <- function(path){
    read.table(paste0(argv$path, "/md5sum.txt"),sep=" ")[,-2]
}
raw_md5 <- do.call(rbind, lapply(word(ln_cmd[,3],1,sep=fixed("/R")), function_read_raw_md5))
colnames(raw_md5) <- c("md5", "sample")
raw_md5 <- raw_md5[!duplicated(raw_md5$md5),]
raw_md5 <- raw_md5[order(raw_md5$md5),]

new_md5 <- read.table(paste0(argv$path, "/md5sum.txt"),sep=" ")[,-2]
colnames(new_md5) <- c("md5", "sample")
new_md5 <- new_md5[order(new_md5$md5),]

# check --
print("md5sum files are equal: ")
print(identical(new_md5$md5, raw_md5$md5))
print("sample id are equal: ")
print(identical(word(new_md5$sample, 2, sep=fixed("/R")), word(raw_md5$sample, 2, sep=fixed("/R"))))