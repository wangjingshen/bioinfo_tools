### istarSpots
该流程用于把 istar 的分群迁移到 spot 上, 并输出文件用于后续分析。


### parameters
--istar_labels        istar labels pickle

--dir                 celescope dir

--spname              spname

--k                   k istar pixels, default = 3, type = int

--distance_thresh     threshold of distance between spots and istar pixels, default=200, type = int

--clip                clip spots from istar, store_true


### note
在图片、像素、空间转录组里, 行 / 列与数学里不一样, 

列（col）= 左右方向 = 水平 = X 轴

行（row）= 上下方向 = 垂直 = Y 轴