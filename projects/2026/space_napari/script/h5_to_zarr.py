# https://spatialdata.scverse.org/en/stable/tutorials/notebooks/notebooks/examples/sdata_from_scratch.html
# celescope_space_outs to spatialdata

import argparse
import scanpy as sc
import pandas as pd
import geopandas as gpd
import spatialdata as sd
from shapely.geometry import Point
from spatialdata.models import Image2DModel, ShapesModel, TableModel
import json
import contextlib
from pathlib import Path
import sys
#import spatialdata_plot  # load spatialdata.pl
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd


def visium_input(space_dir, sample) -> None:
    cmds = [f'cp -r {space_dir}/outs {sample}/',
            f'mv {sample}/spatial/positions_list.csv  {sample}/spatial/tissue_positions_list.csv']
    for i in cmds:
        execute_cmd(i)

def h5_to_zarr(space_dir, sample, outdir) -> None: 
    # read h5
    adata = sc.read_visium(f'{sample}/', count_file="filtered_feature_bc_matrix.h5")
    
    # img
    img = adata.uns['spatial']['library0']['images']['hires']
    img_for_sdata = Image2DModel.parse(data=img, scale_factors=(2, 2, 2), dims=("y", "x", "c"))
    
    with open(f'{sample}/spatial/scalefactors_json.json', 'r') as f:
        scalefactors = json.load(f)
        tissue_hires_scalef = scalefactors['tissue_hires_scalef']
    
    centers = adata.obsm["spatial"] * tissue_hires_scalef  # import, align the scale of spot to scale of H&E hires
    spot_df = pd.concat([adata.obs[["array_col", "array_row"]].reset_index(drop=True), pd.DataFrame(centers, columns=["x", "y"])], axis=1, ignore_index=True,)
    spot_df.columns = ["array_col", "array_row", "spot_center_x", "spot_center_y"]

    fixed_row = spot_df.array_row.iloc[0].item()
    cols = spot_df.query(f"array_row == {fixed_row}").array_col
    min_col, max_col = cols.min().item(), cols.max().item()
    xs = spot_df.query(f"array_row == {fixed_row}").spot_center_x
    min_x, max_x = xs.min().item(), xs.max().item()

    px_per_um = (max_x - min_x) / ((max_col - min_col) / 2) / 100
    # each spot is 55 µm in diameter
    radius = px_per_um * 55 / 2

    # shape
    df = pd.DataFrame([radius] * len(centers), columns=["radius"])
    gdf = gpd.GeoDataFrame(df, geometry=[Point(x, y) for x, y in centers])
    shapes_for_sdata = ShapesModel.parse(gdf)

    # remove remnants of previous way to store spatial data
    with contextlib.suppress(KeyError):
        del adata.uns["spatial"]
    with contextlib.suppress(KeyError):
        del adata.obsm["spatial"]

    adata_for_sdata = TableModel.parse(adata)

    adata_for_sdata.uns["spatialdata_attrs"] = {
        "region": "spots",  
        "region_key": "region",  
        "instance_key": "spot_id",  
    }
    adata.obs["region"] = pd.Categorical(["spots"] * len(adata))
    adata.obs["spot_id"] = shapes_for_sdata.index
    adata.obs[["region", "spot_id"]]

    sdata = sd.SpatialData(
        images={"hires": img_for_sdata},
        shapes={"spots": shapes_for_sdata},
        tables={"adata": adata_for_sdata},
    )

    mkdir(outdir)
    sdata.write(f'{outdir}/{sample}.zarr', overwrite = True)

    #import matplotlib.pyplot as plt
    #fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    #sdata.pl.render_images().pl.render_shapes(color="array_row").pl.show(ax=axs[0], title="Row")
    #sdata.pl.render_images().pl.render_shapes(color="array_col").pl.show(ax=axs[1], title="Col")

    cmd = f'rm -rf {sample}/'
    execute_cmd(cmd)

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--space_dir', help='celescope space dir', required=True)
    parsers.add_argument('--sample', help='sample', required=True)
    parsers.add_argument('--outdir', default = "napari", help='outdir, default:napari')

    args = parsers.parse_args()
    visium_input(args.space_dir, args.sample)
    h5_to_zarr(args.space_dir, args.sample, args.outdir) 

if __name__ == '__main__':
    main()