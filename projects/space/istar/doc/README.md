### istar

该流程用于跑istar的流程


### parameters

--dir                  celescope空转分析目录

--image                H&E染色图

--spname               样本名

--swap_pos             置换 positions_list 文件的最后两列, 即 y = pxl_row_in_fullres， x = pxl_col_in_fullres

--foreground_method    识别前景的方法，默认为 max， 即RGB颜色差异大的区域为前景

--hard_radius          重新设置 radius

--cluster_method       分群方法，默认为 km

--n_cluster            分群个数，默认为 10

--step                 分析步骤，默认为 input,preprocess,get_mask,impute,cluster,downstream,output


### reference

1.https://github.com/daviddaiweizhang/istar