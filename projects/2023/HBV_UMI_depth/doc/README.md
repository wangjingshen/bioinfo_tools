1.概述: 该代码对不同测序深度检测到病毒UMI的数量进行作图。
包含两步: 
1)获取病毒UMI:get_virus_bam.py
2)作图: plot.R

2.需要的包:
python包:pysam,pandas,argparse
R包: tidyverse,argparser

3.运行:
3.1 获取病毒UMI,示例在 test/get_virus_bam
运行代码: 
python path/script/multi.py mapfile   ## 注意：path需要修改为对应的位置

mapfile包含4列, 分别为:
1)富集文库celescope分析路径;
2)UMI的最小支持reads数(从celescope分析路径里的.metrics.json文件中获取);
3)输出路径;
4)结果文件的前缀名;

3.2 作图, 示例在test/test_one_sample及test/test_multi_sample
运行代码:
Rscript path/script/plot.R \
    --virus_id_file ../get_virus_bam/R220124006/R220124006_virus_id.tsv \  ## 注意：path需要修改为对应的位置
    --sample_name R220124006 \
    --total_reads 6147385 \
    --downsample_mode 10

参数说明:
virus_id_file: 病毒UMI文件,从上一步可以获取
sample_name: 样本名
total_reads: 总reads数, 从celescope分析路径里的01.barcode/stat.txt文件中获取
downsample_mode: 下采样的方式, 10(以0.1为步长)或者20(以0.05为步长,默认)
outdir: 输出目录, 默认为当前路径

