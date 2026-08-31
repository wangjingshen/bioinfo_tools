## larry
该脚本用于分析 新格元 v3 的 larry 文库


## run
source activate celescope3.0.0

python ../scripts/pipeline.py \
    --fq_path /SGRNJ06/randd/USER/wangjingshen/script_dev/2026/larry/test/input/ \
    --library_id CUSL260130019 \
    --matrix /SGRNJ06/randd/USER/wangjingshen/script_dev/2026/larry/test/input/larry_p444_EmptyDrops_CR_matrix_10X/larry_p444_EmptyDrops_CR_matrix_10X.tar \
    --spname larry_p444


## 参数
--fq_path         larry 文库 fastq 路径
--library_id      larry 文库 id
--matrix          larry 文库配对的 RNA 矩阵
--spname          larry 样本名


## reference
1.https://github.com/AllonKleinLab/LARRY/blob/master/LARRY_for_10X.ipynb