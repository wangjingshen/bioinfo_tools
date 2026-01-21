# napari doesn’t work on the server, so we’re running it locally.
# recode code here

#import os
#os.environ["PYQTGRAPH_QT_LIB"] = "PyQt5"   # 或 "PySide6"
import pandas as pd
import spatialdata as sd
import spatialdata_plot  # noqa: F401
from napari_spatialdata import Interactive
from spatialdata import polygon_query
import matplotlib.pyplot as plt
# zarr v3 to v2
#import zarr, os
#os.environ["ZARR_V3_EXPERIMENTAL_API"] = "0"
#zarr.config.config["default_format"] = "v2"

# read zarr
visium_sdata = sd.read_zarr("napari/Mint_FFPE_96_96_1119.zarr/")
Interactive(visium_sdata)  # open napari

select_tables = {}
for shape in ["select"]:
    polygon = visium_sdata[shape].geometry.iloc[0]
    select_tables[shape] = polygon_query(visium_sdata, polygon=polygon, target_coordinate_system="global")["adata"]


assert visium_sdata["adata"].obs.index.is_unique     # it's important for the index to be unique
categories = ["spots"] + list(select_tables.keys())
n = len(visium_sdata["adata"])

visium_sdata["adata"].obs["annotation"] = pd.Categorical(["spots" for _ in range(n)], categories=categories)

for shape, subtable in select_tables.items():
    in_shape = subtable.obs.index
    visium_sdata["adata"].obs["annotation"].loc[in_shape] = shape

visium_sdata["adata"].obs["annotation"].value_counts()

plt.figure(figsize=(12, 7))
ax = plt.gca()
visium_sdata.pl.render_images("hires").pl.render_shapes("spots", color="annotation").pl.show(coordinate_systems="global", ax=ax)
plt.savefig('napari/select.png')
plt.savefig('napari/select.pdf', dpi=300, bbox_inches='tight')