# rd script
suppressMessages(suppressWarnings({
    library(tidyverse)
    library(argparser)
    library(ggthemes)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args ----
argv <- arg_parser('')
argv <- add_argument(argv, "--csv", help = "merge csv")
argv <- add_argument(argv, "--spname", help = "spname")
#argv <- add_argument(argv, "--version", help = "v1 or v2")
argv <- add_argument(argv, "--outdir", help = "output dir. Default: outdir")
argv <- parse_args(argv)

csv <- str_split(argv$csv, ",")
if(!is.na(argv$spname)){
    spname <- unlist(str_split(argv$spname, ","))
    print(spname)
}

outdir <- argv$outdir
if(!dir.exists(argv$outdir)){
    dir.create(argv$outdir, recursive = TRUE)
}

sort_merged <- do.call(rbind, lapply(csv, read_csv))
df1<-data.frame(sort_merged)
if(!is.na(argv$spname)){
    df1 <- df1[ df1$Sample_ID %in% spname,]
}

df3 <- data.frame(lapply(df1,function(x) as.numeric(gsub("%", "", x))/100))   # sub to gsub
sample_number <- dim(df1)[1]
#png(file = 'RawReads.png')
p1 <- ggplot(df1,aes(Sample_ID,Raw_Reads,fill=Raw_Reads)) +
    geom_bar(position = "dodge",stat = "identity",width = 0.5)+ ##ggplot2的几何图像层（至少要到这一层）
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5, vjust=1, size=16)) + 
    xlab("SampleName") + 
    ylab('Datasize')+
    geom_text(aes(y= Raw_Reads*0.99,label=Raw_Reads),position = position_dodge(0.9), vjust = -0.5 )+  #添加标签
    ggtitle(label = "Rawreads summary")
##mapping_plot

#堆积柱形图
sample_name <- rep(as.vector(df1$Sample_ID),3)

map_group <- c(rep("Mapped_Reads_Assigned_To_Exonic_Regions",sample_number),  # add Mapped_
               rep("Mapped_Reads_Assigned_To_Intronic_Regions",sample_number),
               rep("Mapped_Reads_Assigned_To_Intergenic_Regions",sample_number))

if(sum(grepl("Mapped_Reads_Assigned_To_Exonic_Regions",colnames(df3)))){
    value1<-c(as.vector(df3$Mapped_Reads_Assigned_To_Exonic_Regions),  # add Mapped_
              as.vector(df3$Mapped_Reads_Assigned_To_Intronic_Regions),
              as.vector(df3$Mapped_Reads_Assigned_To_Intergenic_Regions))
}else{
    value1<-c(as.vector(df3$Reads_Assigned_To_Exonic_Regions),  # add Mapped_
              as.vector(df3$Reads_Assigned_To_Intronic_Regions),
              as.vector(df3$Reads_Assigned_To_Intergenic_Regions))   
}

df2 <- data.frame(sample_name,
                  map_group,
                  value1)
#png(file = 'map_region.png')
ggplot(df2, aes(x=sample_name, y=value1, fill=map_group))+
    geom_bar(stat="identity",width = 0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),plot.title = element_text(hjust = 0.5, vjust=1, size=16)) + 
    ylab("Mapping rate") +
    scale_y_continuous(labels = scales::percent)+
    scale_fill_brewer(palette = 'Set2') + #设置颜色板
    geom_text(aes(label=scales::percent(value1,accuracy = 0.01)),position=position_stack(vjust=0.5),size=3.5)+#添加标签
    ggtitle(label = "Mapping Summary")
ggsave(str_glue("{outdir}/mapping_summary.png"), width =8, height=6)


###重要测序指标一览表
#sample_number<-4
sample_name <- rep(as.vector(df1$Sample_ID),4)
quota_group <- c(rep("ValidReads",sample_number),
                 rep("UniquelyMappedReads",sample_number),
                 rep("Fraction_Reads_in_Cells",sample_number),
                 rep("Saturation",sample_number))

value2<-c(as.vector(df3$Valid_Reads),
          as.vector(df3$Reads_Mapped_Uniquely_To_Transcriptome),
          as.vector(df3$Fraction_Reads_in_Cells),
          as.vector(df3$Saturation))

df4<-data.frame(sample_name,quota_group,value2)

p3 <- ggplot(df4,aes(x=quota_group,y=value2,fill=sample_name))+
  geom_bar(position = "dodge",stat = "identity",width = 0.9)+ 
  theme(axis.text.x = element_text( hjust = 0.5, vjust = 0.5),plot.title = element_text(hjust = 0.5, vjust=1, size=16)) +
  xlab('')+
  scale_y_continuous(labels = scales::percent)+
  ylab("percentage of quota") +
  geom_text(aes(y= value2-0.05,label=scales::percent(value2,accuracy = 0.01)),position = position_dodge(0.9), vjust = 0.5 )+  #添加标签
  ggtitle(label = "Important quota")+
  coord_flip() #旋转坐标轴


##cell summary
sample_name <- rep(as.vector(df1$Sample_ID),5)
cell_summary_group <- c(rep("EstimatedNumberofCells",sample_number),rep("MedianUMIperCell",sample_number),
                        rep("TotalGenes",sample_number),rep("MedianGenesperCell",sample_number),rep("MeanReadsperCell",sample_number))
cell_summary_group
value3 <- c(as.vector(df3$Estimated_Number_of_Cells),
            as.vector(df3$Median_UMI_per_Cell),
            as.vector(df3$Total_Genes),
            as.vector(df3$Median_Genes_per_Cell),
            as.vector(df3$Mean_Reads_per_Cell))*100

df5 <- data.frame(sample_name,
                  cell_summary_group,
                  value3)

ggplot(df5,aes(x=cell_summary_group,y=value3,fill=sample_name))+
    geom_bar(position = "dodge",stat = "identity",width = 0.9)+ 
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), plot.title = element_text(hjust = 0.5, vjust=1, size=16)) +
    xlab('')+
    ylab("number") +
    ggtitle(label = "Cell summary")+
    geom_text(aes(y= value3+0.05,label=value3),position = position_dodge(0.9), vjust = -0.5, size=3 )  #添加标签
ggsave(str_glue("{outdir}/important_quota.png"), width =8, height=6)

