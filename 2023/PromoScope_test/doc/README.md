1.糖基化项目分析需求:
1) 糖基化丰度展示;
2) 高、低糖基化丰度细胞的比较分析;
注:后续根据客户的新需求进行更新;


2.分析脚本:
/SGRNJ06/randd/USER/wangjingshen/script/sweetseq/script/analysis.R

1)参数说明:
--rds                注释好的rds, 这里输入的是单一样本的rds, 不同样本的糖基化丰度不可比较;
--sweet_tsne_tag     糖基化分析目录中的丰度文件 (*/04.count_tag/*_tsne_tag.tsv, *需要改为具体路径和文件名);
--species            物种, 人或小鼠, 默认为人;
--filter_cluster     需要过滤掉的细胞簇名, 默认去除doublet;
--max_cutoff         糖基化丰度的可视化上限, 默认为10000;
--analysis_cluster   用于比较分析的细胞簇名, 多个细胞簇名以逗号分隔, 默认对细胞数最多的簇进行分析;
--percent            糖基化丰度分组的比例, 默认为0.2, 即取丰度最高的20%的细胞作为高丰度组, 取丰度最低的20%的细胞作为低丰度组;
--outdir             输出目录, 默认为outdir;
--width_go_high      高丰度组go富集分析的结果图的宽, 默认为10;
--height_go_high     高丰度组go富集分析的结果图的高, 默认为8;
--width_go_low       低丰度组go富集分析的结果图的宽, 默认为10;
--height_go_low      低丰度组go富集分析的结果图的高, 默认为8;
注: 富集分析结果图可使用输出目录中的xxx_go_res.rds进行调整

2)运行示例:
/SGRNJ06/randd/USER/wangjingshen/script/sweetseq/test/test/run.sh
请注意:测试数据为客户样本, 请勿外传;
