import cv2
import json
import numpy as np
import pandas as pd

def get_mask(dir, spname):
    img = cv2.imread(f"{dir}/outs/spatial/tissue_hires_image.png")
    scale = json.load(open(f"{dir}/outs/spatial/scalefactors_json.json"))["tissue_hires_scalef"]

    spot_px = pd.read_csv(f"{dir}/outs/spatial/positions_list.csv", header=None, index_col=0)
    spot_px = spot_px[spot_px[1] == 1][[5, 4]] * scale  # in_tissue=1

    print(spot_px)

    mask = np.zeros(img.shape[:2], dtype=np.uint8)
    for x, y in spot_px.values:
        cv2.circle(mask, (int(x), int(y)), int(scale * 0.5), 255, -1)

    bg_img = cv2.bitwise_and(img, img, mask=~mask)
    cv2.imwrite(f"{spname}/mask-small.png", bg_img)