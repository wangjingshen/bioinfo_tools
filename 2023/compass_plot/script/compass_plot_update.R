suppressWarnings(suppressMessages({
    library(argparser)
    library(data.table)
    library(tidyverse)
    library(ggrepel)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--wilcox_res", help = "wilcox result from running compass_plot.py")
argv <- add_argument(argv, "--subsystem", help = "subsystem, need to be consistent with the subsystem parameter used for running compass_plot.py")
argv <- add_argument(argv, "--xlab_name", help = "xlab name of compass_plot.pdf from compass_plot.py")
argv <- add_argument(argv, "--annotation_text_size", help = "annotation text size, default: 3")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- add_argument(argv, "--name", help = "name, need to be consistent with the name parameter used for running compass_plot.py")
argv <- parse_args(argv)

if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir) 
}
argv$annotation_text_size <- ifelse(is.na(argv$annotation_text_size), 3, as.numeric(argv$annotation_text_size))

# read table --
wilcox_res <- fread(argv$wilcox_res)
plot_data <- wilcox_res %>%
    filter(subsystem == argv$subsystem)
p1 <- plot_data %>%
    ggplot(aes(x = cohens_d, y = -log10(adjusted_pval))) +
    geom_point() +
    xlim(c(-max(abs(plot_data$cohens_d)), max(abs(plot_data$cohens_d)))) +
    xlab(paste0("Cohen's d (", argv$xlab_name, ")"))+
    ylab('-log10(adjusted_pval)') +
    ggtitle(argv$subsystem) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(argv$outdir, "/", argv$name, "_update.pdf"), plot=p1, height=5, width=6)
ggsave(paste0(argv$outdir, "/", argv$name, "_update.png"), plot=p1, height=5, width=6)

p2 <- plot_data %>%
    ggplot(aes(x = cohens_d, y = -log10(adjusted_pval))) +
    geom_text_repel(aes(label = reaction_name), size = argv$annotation_text_size, max.overlaps =30) +
    geom_point() +
    xlim(c(-max(abs(plot_data$cohens_d)), max(abs(plot_data$cohens_d)))) +
    xlab(paste0("Cohen's d (", argv$xlab_name, ")"))+
    ylab('-log10(adjusted_pval)') +
    ggtitle(argv$subsystem) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(argv$outdir, "/", argv$name, "_anno_update.pdf"), plot=p2, height=10, width=12)
ggsave(paste0(argv$outdir, "/", argv$name, "_anno_update.png"), plot=p2, height=10, width=12)


print("Compass plot update done.")