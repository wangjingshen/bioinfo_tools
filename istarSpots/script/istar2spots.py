import sys
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from pathlib import Path

ROOT = Path(__file__).resolve().parent
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, timer, load_pickle


def map_spot_to_istar(spot_x, spot_y, kdtree, istar_clusters, k = 5, distance_thresh = 200, clip = False):
    '''
    most cluster in k istar pixels

    clip TRUE: The result is similar to istar
    clip FALSE: The result is similar to space (recommend)
    '''
    if not (0 <= spot_x <= 1008 and 0 <= spot_y <= 1008):   # 0 <= spot_x < 1008 and 0 <= spot_y < 1008
        if(clip == True):
            logger.info(f"spot {spot_x},{spot_y} beyond the istar area, clipped.")
            return -1
        else:
            logger.info(f"spot {spot_x},{spot_y} beyond the istar area, retained.")

    distances, indices = kdtree.query([spot_x, spot_y], k=k)
    valid_idx = indices[distances < distance_thresh]
    if len(valid_idx) == 0:
        logger.info(f"spot {spot_x},{spot_y} invalid, consider adjusting distance_thresh.")
        return -1
    nearest_clusters = istar_clusters[valid_idx]
    return np.bincount(nearest_clusters).argmax()


def istar2spots(istar_lables, spots_pos, outdir, k, distance_thresh, clip):
    istar_cluster = load_pickle(istar_lables)  # (1008, 1008)
    spot_coords = pd.read_csv(spots_pos, header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"])  # "pxl_col", "pxl_row"
    spot_coords = spot_coords[spot_coords["in_tissue"] == 1].reset_index(drop=True)

    # scale 1008*1008
    spot_coords["pxl_col_abs"] = spot_coords["pxl_col"] - spot_coords["pxl_col"].min()
    spot_coords["pxl_row_abs"] = spot_coords["pxl_row"] - spot_coords["pxl_row"].min()
    scale_col = 1008 / spot_coords["pxl_col_abs"].max()
    scale_row = 1008 / spot_coords["pxl_row_abs"].max()
    spot_coords["x"] = (spot_coords["pxl_col_abs"] * scale_col).astype(int)
    spot_coords["y"] = (spot_coords["pxl_row_abs"] * scale_row).astype(int)

    y, x = np.where(istar_cluster != -1)
    istar_coords = np.column_stack([x, y])
    istar_labels = istar_cluster[y, x]
    kdtree = KDTree(istar_coords)

    spot_coords["cluster"] = spot_coords.apply(
        lambda r: map_spot_to_istar(r["x"], r["y"], kdtree, istar_labels, k=k, distance_thresh=distance_thresh, clip=clip), axis=1
    )

    valid = spot_coords[spot_coords["cluster"] != -1]

    plt.figure(figsize=(9, 8))
    #plt.imshow(istar_cluster, cmap="tab20", alpha=0.2, origin="upper")  # show istar raw plot 
    plt.scatter(
        valid["x"], valid["y"],
        c=valid["cluster"],
        cmap="tab10",
        s=12,
        edgecolor="black",
        lw=0.2
    )
    plt.gca().invert_yaxis()  # 
    plt.colorbar(shrink=0.7)
    plt.axis("off")
    plt.title("iSTAR", fontsize=13)
    plt.savefig(f"{outdir}/istar2spots.png", dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close()

    result_df = spot_coords[["barcode", "cluster"]].copy()
    result_df.columns = ["barcode", "istar_cluster"]
    result_df.to_csv(f'{outdir}/istar2spots_df.csv', index=False)
    #spot_coords.to_csv(f'{spname}/istar2spots_df_total.csv', index=False)

    logger.info(f"valid Spot: {len(valid)} / {len(spot_coords)}")
    logger.info(valid["cluster"].value_counts().sort_index())

    #with open(f'{outdir}/summary.txt', "w") as f:
    #    f.write(f"valid Spot: {len(valid)} / {len(spot_coords)}\n")
    #    f.write("\nistar cluster:\n")
    #    f.write(valid["cluster"].value_counts().sort_index().to_string())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--istar_labels', help='istar labels pickle', required=True)
    parser.add_argument('--space_pos', help='space positions', required=True)
    parser.add_argument('--outdir', help='outdir', required=True)
    parser.add_argument('--k', default = 3, type = int, help='k istar pixels')
    parser.add_argument('--distance_thresh', default=200, type = int, help='thresh of distance between spots and istar pixels')
    parser.add_argument('--clip', action='store_true', help='clip spots from istar') 
    args = parser.parse_args()

    istar2spots(args.istar_labels, args.space_pos, args.outdir, args.k, args.distance_thresh, args.clip)


if __name__ == '__main__':
    main()