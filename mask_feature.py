import os
from pathlib import Path
import rasterio
from rasterio.io import MemoryFile 
from rasterio import features
from rasterio.enums import Resampling
from rasterio import Affine
import arcpy
import json
import numpy as np
from shapely.geometry import shape
import time


#user_dir = os.path.expanduser('~')

#script_path = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 'anaconda3', 'Scripts')
#gdal_data = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 'anaconda3', 'envs', 'LiDOG', 'Library', 'share', 'gdal')
#proj_lib =Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 'anaconda3', 'envs', 'LiDOG', 'Library', 'share')

#if script_path.name not in os.environ["PATH"]:
#    os.environ["PATH"] += os.pathsep + str(script_path)

#os.environ["GDAL_DATA"] = str(gdal_data)
#os.environ["PROJ_LIB"] = str(proj_lib)


r_path = r'X:\iocm_deliverables\ocs\FL1609\dem\FL1609-TB-N_mllw_dem.tif'
r_coverage_geojson = r'X:\iocm_deliverables\ocs\FL1609\dem\FL1609-TB-N_mllw_dem_COVERAGE.geojson'
r_coverage_shp = r'X:\iocm_deliverables\ocs\FL1609\dem\FL1609-TB-N_mllw_dem_COVERAGE.shp'

print('determining coverage of DEM...')

crs = {
        "type": "name",
        "properties": {
          "name": "epsg:26916"
        }
      }

tic = time.time()

with rasterio.open(r_path) as r:
    t=r.meta['transform']
    print(r.crs.to_epsg())
    data = r.read_masks(1, out_shape=(r.height // 10, r.width // 10), resampling=Resampling.average)


transform = Affine(t.a * 10, t.b, t.c, t.d, t.e * 10, t.f)
data = np.ma.array(data, mask=(data == 0))
mask = features.shapes(data, mask=None, transform=transform)

print('creating in_memory src dem coverage...')
coverage_shp = 'in_memory\cov'
sr = arcpy.SpatialReference(r.crs.to_epsg())
arcpy.CreateFeatureclass_management('in_memory', 'cov', spatial_reference=sr)
cursor = arcpy.da.InsertCursor(coverage_shp, ['SHAPE@WKT'])

print('adding polygons to coverage shp...')
for poly, value in mask:
    poly_wkt = shape(poly).wkt
    cursor.insertRow([poly_wkt])

del cursor

print(time.time() - tic)
