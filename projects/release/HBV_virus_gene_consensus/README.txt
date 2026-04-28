1.概述: 该代码用于对判定属于病毒UMI的reads进行consensus.

2.需要的包:
celescope1.5.1b0
pysam,pandas,argparse,biopython,xopen

3.运行:
示例在test/

运行代码: 
python path/multi.py mapfile   ## 注意：path需要修改为对应的位置


mafile包含7列,分别为:
1)code_path: analysis.py所在路径;
2)fj_path: 富集文库celescope分析路径;
3)max_n_reads: UMI的reads数阈值,小于等于阈值的先进行多序列比对再进行consensus,大于阈值的使用celescope内置的consensus;
4)small_threshold: 先进行多序列比对再进行consensus方法的最常见碱基的频率阈值,默认为0.5,超过该阈值设置为该碱基,否则设置为'N';
5)large_threshold: 使用celecope内置的consensus方法的最常见碱基的频率阈值,默认为0.5,超过该阈值设置为该碱基,否则设置为'N';
6)min_consensus_read: 最常见碱基的最少出现reads数,默认为1,即最少出现1次;
7)outdir: 输出目录.