## 前言
此脚本用于 celescope FFPE 识别 circRNA，分析基于 CIRCexplorer2 [1]。

GEXSCOPE FFPE 单细胞核转录组（新格元）采用随机引物原位反转录，不依赖 polyA 捕获，属于无偏全转录组策略。和基于 oligo-dT 的 polyA 单细胞转录组不同，理论上可以捕获无 polyA 尾的环状 RNA 片段，具备检测 circRNA 的基础可能性。但该体系存在两个"短板"，会直接限制 BSJ（反向剪接断点）的检出：1）核内 circRNA 本底丰度低，绝大多数外显子来源 circRNA 成熟后会转运至细胞质，细胞核内仅少量内含子 ciRNA 及部分滞留核内的外显子型 circRNA，因此天然丢失大量胞质环状 RNA；2）FFPE 交联导致 RNA 断裂，容易破坏 BSJ 断点：福尔马林固定带来的交联以及后续脱交联加热处理会打断 RNA。即便随机引物扩增到 circRNA 的部分序列，一旦反向剪接断点被打断，就无法生成跨越 BSJ 的嵌合 reads，STAR 和 CIRCexplorer2 也就不能识别出 circRNA。蜡块储存时间越长、RNA 的 DV200 （长度＞200 nt 的 RNA 片段，占全部 RNA 片段的百分比）越低，BSJ 的丢失情况会越明显。

## 分析步骤
step1.使用 STAR 识别嵌合转录本, 参数参照已发表文献[2]
step2.使用 CIRCexplorer2 鉴定环状 RNA
step3.将 CIRCexplorer2 结果分配到 barcode 和 UMI

## 运行
```
source activate celescope3.0.0

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/rna/circexplorer2/scripts/get_celescope_sjm.py \
    --mapfile mapfile \
    --product ffpe \
    --genome /SGRNJ06/randd/public/genome/rna/mmu/mmu_ensembl_110_nofilter \
```
投递 sjm 命令，等运行完成再运行下面的代码

```
source activate celescope3.0.0

python /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/rna/circexplorer2/scripts/pipeline.py \
    --celescope_dir /SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/rna/circexplorer2/test/mouse-testicle/ \
    --name mouse-testicle \
    --refFlat /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/circexplorer2/genome/Mus_musculus.GRCm39.110.refFlat_11col \
    --reference_fa /SGRNJ06/randd/USER/wangjingshen/rd_project/2026/circexplorer2/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --match_window 10
```

## 参数
get_celescope_sjm.py 大部分参数和 celescope 一致，额外的参数主要和 STAR 嵌合检测相关，目前只添加了 outFilterMismatchNmax， 后续根据测试结果考虑纳入其他相关参数。

pipeline.py 参数

--celescope_dir  celescope 分析目录
--name           样本名
--refFlat        参考基因组 refFlat 文件，需要 11 列
--reference_fa   参考基因组 fasta 文件
--match_window   circexplorer2注释环状RNA会修正坐标，这里基于 10 window 将 注释信息 和 barcode、UMI 进行匹配

### STAR 嵌合检测参数说明（circRNA BSJ 识别）
--chimSegmentMin：嵌合片段最小比对长度，默认值为 0，设为 0 时 STAR 不会输出任何嵌合比对结果；本分析参考文献设置为 20。嵌合 reads 会拆分为 A、B 两段分别比对基因组，只有 A、B 两段各自的比对长度均不小于设定值，该 read 才会被判定为嵌合 read。

--chimScoreMin：嵌合两段比对总分的最低阈值，总分低于阈值的嵌合 read 会被 STAR 丢弃，不会写入 Chimeric.out.junction，默认值 0；本分析参考文献设置为 1。嵌合 read 拆分为 A、B 两段分别比对（BSJ 属于典型场景），嵌合总得分等于 A 段比对得分加 B 段比对得分，单段比对得分 = 匹配碱基数 − 错配罚分 − gap 罚分。

--alignIntronMax：线性剪接最大允许内含子长度，默认 0，由 STAR 自动计算，计算公式：maxIntron = (2^winBinNbits) × winAnchorDistNbins。默认 winBinNbits=14，winAnchorDistNbins=9，自动计算得到最大内含子约 147456 bp。该参数仅影响普通 mRNA 线性剪接比对，不会改变 circRNA 的 BSJ 检出数量，本分析与文献保持一致，设置为 1000000。

--outFilterMismatchNmax：单条 read 比对允许的最大错配总数，默认 10，错配数量超过阈值则丢弃这条比对。文献采用数值 4，相比默认值更为严格。(FFPE 样本存在 RNA 降解与碱基化学修饰，碱基错误率更高，很多真实 BSJ 对应的 reads 错配落在 4~10 之间，会被该参数过滤)

--alignTranscriptsPerReadNmax：单条 read 最多检索的比对构型数量，默认 10000；达到上限就停止继续查找更多比对方案，防止提前截断而丢失 BSJ 候选比对。

--outFilterMultimapNmax：一条 read 允许比对到基因组的最大位点数。若候选定位位点超过阈值，整条 read 直接丢弃，线性比对与嵌合检测都会受到影响，默认 10。该参数针对完整 read 的基因组多定位，属于前置过滤：如果 read 能比对到超过 2 个基因组位置，会直接丢弃，没有机会进入嵌合检测流程，间接影响 BSJ 检出。

--chimMultimapNmax：嵌合模式下允许的嵌合片段多比对结果数量，默认 0。设置为 0 时启用旧版嵌合算法，仅保留两段均为唯一比对的 read；A 段或 B 段任意一个嵌合片段存在多 mapping，该嵌合 junction 直接丢弃。

--chimOutType Junctions：输出 Chimeric.out.junction 文本文件，保存全部嵌合断点信息，该文件作为 CIRCexplorer2 的输入，用于识别 BSJ。

--chimJunctionOverhangMin：嵌合断点两侧远端锚定 overhang 的最小长度。跨断点的两段比对，远离断点一侧的有效比对长度均不能低于该阈值。嵌合 read 在 BSJ 断点处拆分为 A、B 两段，靠近断点一侧为 breakpoint，远离断点的末端区域即为 overhang；设置为 20 时，A、B 两段的 overhang 都≥20bp，这条嵌合 junction 才会被保留。


### STAR 嵌合检测逻辑(基于参考文章分析参数)
STAR 首先搜索 read 的候选比对构型，由 alignTranscriptsPerReadNmax 控制最多查找的比对方案数量。随后进行两轮前置过滤：如果一条 read 在基因组上可比对到超过 2 个不同位置，会直接丢弃；保留的 read 再检查总错配数，错配超过 4 个也会被丢弃。

通过前置过滤的 read 会进入嵌合检测流程，STAR 尝试将 read 拆分为两段做嵌合比对。之后依次进行四层嵌合条件校验：两段嵌合片段各自的比对长度不少于 20bp，断点两侧远端锚定区域长度均不少于 20bp，两段比对分数相加不低于 1，并且拆分后的两个嵌合片段都必须是基因组唯一比对。只有全部条件都满足，这条嵌合断点才会输出到 Chimeric.out.junction 文件。

具体而言
Step1：比对候选搜索阶段
参数：--alignTranscriptsPerReadNmax 10000
STAR 最多搜索 10000 种 read 的比对构型，一旦达到上限就停止检索，不再查找更多候选比对方案。该设置用于避免比对搜索提前终止，防止遗漏 BSJ 嵌合的候选比对。本步骤不会丢弃 read，仅限制候选比对的搜索数量。

Step2：整条 read 多定位前置过滤（线性层面）
参数：--outFilterMultimapNmax 2
若一条 read 可比对到基因组上超过 2 个不同位置，该 read 会被直接丢弃，不再参与后续比对和嵌合检测。只有基因组定位位点数量小于或等于 2 的 read，才会进入下一阶段处理。

Step3：整条 read 总错配过滤
参数：--outFilterMismatchNmax 4
评估该 read 最优比对结果的错配数量，当错配碱基总数大于 4 时，丢弃这条 read。错配数不超过 4 的 read 可继续后续流程；插入缺失（gap/indel）不计入错配统计，但会影响比对打分。

Step4：线性剪接内含子长度限制
参数：--alignIntronMax 1000000
该参数仅针对正常 mRNA 的线性剪接，限定线性内含子最大长度为 1 Mb。环状 RNA 反向剪接产生的基因组间隔不受该参数约束，不会对 BSJ 断点进行筛选。

Step5：嵌合检测与 chim 系列参数逐级过滤
只有通过前面所有过滤条件的 read，STAR 才会尝试将 read 拆分为 A、B 两段进行嵌合比对，并依次执行四层校验：

Step5-1 嵌合片段最小比对长度
参数：--chimSegmentMin 20
拆分得到的 A、B 两个嵌合片段，各自有效比对长度均不能少于 20 bp；任意一段比对长度不足 20 bp，即丢弃该嵌合 read。

Step5-2 嵌合断点远端锚定长度
参数：--chimJunctionOverhangMin 20
A、B 两段远离断点一侧的锚定区域（overhang），都需要至少 20 bp。任意一端锚定片段长度不足，则丢弃这条嵌合 junction。

Step5-3 嵌合比对总分阈值

参数：--chimScoreMin 1
A 段比对得分与 B 段比对得分之和不能低于 1，总分小于 1 则丢弃。比对得分计算规则：匹配碱基数 − 错配罚分 − gap 罚分。

Step5-4 嵌合片段多比对控制

参数：--chimMultimapNmax 0
参数设为 0 时启用旧版嵌合检测算法，要求 A 片段和 B 片段均为基因组唯一比对。只要 A、B 任意一个嵌合片段存在多 mapping，该嵌合 junction 直接丢弃。

Step6：结果输出

参数：--chimOutType Junctions
全部通过上述所有过滤条件的嵌合断点信息，最终输出至 Chimeric.out.junction 文件。


### 参数优化
基于已发表文献参数进行测试未识别到带 barcode 和 umi 的嵌合转录本，因为 FFPE 会限制环状 RNA 的检出，考虑放宽参数限制。

outFilterMismatchNmax  从 4 放宽到 10，依旧未检出

chimScoreMin           从 1 放宽到 0，依旧未检出

chimMultimapNmax       从 0 放宽到 1，可以检出

chimMultimapNmax 管控被 STAR 拆分后的两段嵌合片段的多重比对，与全局参数 outFilterMultimapNmax 属于两套独立过滤体系。

当 chimMultimapNmax 0 时，启用 STAR 旧版 BSJ 检测逻辑：嵌合读段拆分得到的两个片段，都要求为唯一比对（每个片段仅能匹配基因组上 1 个位置）。一旦 A 片段或者 B 片段任意一段存在多比对情况，这条嵌合 read 就会被直接丢弃。
优点：原始检出结果背景噪音低，假阳性少；
缺点：FFPE 样本 RNA 存在降解、碱基损伤，很多片段会落入多比对区域，BSJ 检出灵敏度偏低，容易出现检出数量少甚至无结果。

当 chimMultimapNmax 设为 1 时，允许嵌合片段最多 1 个多重比对位点：两段嵌合片段中，最多其中一段可以存在多比对。放宽了嵌合片段的比对唯一性限制，能保留更多 BSJ 候选 read。
优点：提升 BSJ 检出灵敏度，适配 FFPE 降解样本；
缺点：候选集中混入更多噪音，假阳性风险会上升。



## reference
1.https://circexplorer2.readthedocs.io/en/latest/

2.Expression profiling and in situ screening of circular RNAs in human tissues. Sci Rep. 2018 Nov 16;8(1):16953.