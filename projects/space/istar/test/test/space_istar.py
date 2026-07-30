import cv2
import numpy as np

# 读取图片（换成你的第一张图路径）
img = cv2.imread("/SGRNJ06/randd/USER/wangjingshen/rd_project/2026/space/r3/mouse_intestine/Mint_FFPE_96_96_1119/outs/space_pathseq/space_seurat_clusters.png")
hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

# 创建空白画布
result = np.zeros_like(img)

# ====================== 颜色区间（匹配你的图）======================
colors = [
    ([0, 40, 40], [10, 255, 255], (0, 0, 255)),       # 红
    ([11, 40, 40], [25, 255, 255], (0, 127, 255)),   # 橙
    ([26, 40, 40], [35, 255, 255], (0, 255, 255)),    # 黄
    ([36, 40, 40], [70, 255, 255], (0, 255, 0)),      # 绿
    ([71, 40, 40], [100, 255, 255], (255, 127, 0)),  # 青
    ([101, 40, 40], [130, 255, 255], (255, 0, 0)),    # 蓝
    ([131, 40, 40], [165, 255, 255], (127, 0, 127)), # 紫
]

# ====================== 核心：膨胀连接 ======================
kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (12, 12))

for hsv_low, hsv_high, bgr in colors:
    # 提取颜色
    mask = cv2.inRange(hsv, np.array(hsv_low), np.array(hsv_high))
    
    # 多次膨胀 → 把点连成片
    mask = cv2.dilate(mask, kernel, iterations=6)
    
    # 上色
    result[mask > 0] = bgr

# ====================== 输出最终图 ======================
cv2.imwrite("output.png", result)
cv2.imshow("result", result)
cv2.waitKey(0)