该文档用于说明运行 celescope pathseq。

conda环境： celescope2.2.0

## step1.生成 sjm 脚本

脚本示例：
```
## 示例路径：/SGRNJ06/randd/USER/wangjingshen/project/P24120305_BDKQ/B1/r1_110/run.sh
source activate celescope2.2.0

multi_pathseq \
 --mapfile mapfile \
 --genomeDir /SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/Genome/celescope_v2/human_GRCh38_110 \
 --filter_bwa_image /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.fa.img \
 --kmer_file /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.bfi \
 --microbe_bwa_image /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.fa.img \
 --microbe_dict /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.dict \
 --microbe_taxonomy_file /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_taxonomy.db \
 --downsample_reads 10000000 \
 --mod sjm

# 更新 sjm 脚本设置更高的 unlimit 参数，一定要加上下面这一行
python /SGRNJ06/randd/USER/wangjingshen/script/pathseq_sjm_update/script/pathseq_sjm_update.py --q product.q
```

参考基因组说明：

celescope pathseq 的参考基因组分为两部分，宿主基因组（genomeDir，filter_bwa_image，kmer_file）和细菌基因组（microbe_bwa_image，microbe_taxonomy_file，microbe_taxonomy_file），分析中只需要修改宿主基因组。

#### genomeDir
genomeDir 按运营需求，需要和转录组用的参考基因组版本匹配，运营会在邮件中标注版本，对应位置如下。

新格元-v1-92		/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92
                    /SGRNJ/Public/Database/genome/mus_musculus/ensembl_92

新格元-v1-99	    /SGRNJ06/randd/public/genome/rna/hs/hs_ensembl_99
                    /SGRNJ06/randd/public/genome/rna/mmu/mmu_ensembl_99

新格元-v2-99	    /SGRNJ06/randd/public/genome/rna/celescope_v2/hs
                    /SGRNJ06/randd/public/genome/rna/celescope_v2/mmu

新格元-v2-110		/SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/Genome/celescope_v2/human_GRCh38_110
                    /SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/Genome/celescope_v2/mmu_ensembl_110

#### filter_bwa_image，kmer_file
这两个参数是 pathseq 过滤宿主 reads 使用的宿主文件，分为人和小鼠（其余物种需要联系王靖深生成）。

人：  --filter_bwa_image /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.fa.img \
      --kmer_file /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.bfi \

鼠：  --filter_bwa_image /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm10/mm10.fasta.img \
      --kmer_file /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm10/mm10.bfi \


## step2.投递任务
sjm/sjm.job


## step3.重新投递报错的任务(可选)
部分样本pathseq那一步会报错，需要重新投递试一下，命令如下：

cp sjm/sjm.job.status sjm/sjm.job
sjm sjm/sjm.job

注：重新投递仍报错，和王靖深联系进行排查。


## step4.研发这的初步分析
等16S 富集文库跑完，可以跑一些初步分析，包括 转录组的自动注释，细菌检出细胞数，top10检出菌的映射等，脚本如下：

```
## 示例路径：/SGRNJ06/randd/USER/wangjingshen/project/P24121206_ZRYH/B1/r1_110/analysis/run.sh
source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/script/sc16S_analysis_pipeline/script/pipeline.py \
    --matrix_10X /SGRNJ06/randd/USER/wangjingshen/project/P24121206_ZRYH/B1/r1_110/zl/DG/outs/filtered/,/SGRNJ06/randd/USER/wangjingshen/project/P24121206_ZRYH/B1/r1_110/zl/GG/outs/filtered/ \
    --species mouse \
    --rna_spname DG_RNA,GG_RNA \
    --gname DG_RNA,GG_RNA \
    --fj_path /SGRNJ06/randd/USER/wangjingshen/project/P24121206_ZRYH/B1/r1_110/ \
    --fj_spname DG_Target,GG_Target \
```

主要参数：
--matrix_10X: 转录组的mtx路径
--species：human or mouse， 用于自动注释
--rna_spname: 转录组的样本名
--gname： 转录组的组名
--fj_path： 16S富集分析路径
--fj_spname：  16S富集文库名

