## update log

2026.02.28  添加参数mask_method，用于解决mask识别错误问题，有的H&E图会呈现背景差异大的情况；

2026.03.26  去掉参数mapfile，改为具体参数，后面考虑增加 multi 以支持多样本运行；
            istar中impute 可能出现空的patch，修复 1）使用0补齐；2）使用 hard_radius（调小）（需要评估一下效果）；
            重命名参数 mask_method 为 foreground_method，避免误解；
