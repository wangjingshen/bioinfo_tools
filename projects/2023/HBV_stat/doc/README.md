1.get_prefix.R
该脚本获取转录组细胞的前缀名，用于与富集文库进行匹配;
注：特别针对走lims分析的rds，lims分析会对细胞进行质量过滤.


2. HBV_stat.R
该脚本用于统计每个细胞类型中HBV的检出率；

参数:
--rds, rds路径, 可以 1)输入lims审核的rds(推荐),2)邮件审核的rds(弃用,仅作为支持此前结果的复现),3)未注释的rds;三者的区别见注1
--version: rds的版本,v1或者v2(默认); v1是邮件审核的rds; v2是lims审核的rds,
--anno_file: 注释文件,该参数针对未注释的rds,输入文件参照auto_assign的结果,需要两列,一列为 "cluster", 从1开始; 第二列为 "cell_type",即每个cluster对应的细胞名称;示例见注2
--cele_version:  celescope版本,不同版本的富集文库目录结构、命名有所差别,目前主要针对1.5.1(默认)以及1.13.0;
--mod:  富集文库目录模式, same或者no(默认); same是指所有富集文库都在同一目录下, no是指所有富集文库并非在同一目录下;
--fj_path:  富集文库路径, 如果mod参数设置为same,那么此处输入一个路径即可; 如果mod参数设置为no,则需要输入所有的富集文库路径,使用逗号进行分隔;
--fj_name:  富集文库名称, 使用逗号进行分隔;
--prefix:  转录组细胞名前缀,一般是转录组样本名,不清楚的话可以使用get_prefix.R获取前缀名;
--subset:  该参数用于从整合注释的rds中选取部分样本进行分析, 输入要选取的样本名,使用逗号分隔;
--outdir:  输出目录,默认为 HBV_stat/ .

注1: lims审核的rds和邮件审核的rds的主要区别在于:
1) lims审核的rds进行了质量过滤,因此富集文库也需要进行过滤;
2) 另外,lims审核的rds的细胞注释变量是cluster,样本名称是sample;
        而邮件审核的rds的细胞注释变量是cell_type,样本名称是orig.ident.

注2: anno_file示例:
cluster cell_type
1   T
2   B




版本说明:
V0: 初始版本;
V1: 取消设置默认阳性细胞在上;
V2: 输出文件名设置为prefix
V3(当前使用):自动识别多版本的celescope, 优化一些写法