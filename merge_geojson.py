from pathlib import Path
from functools import partial
import json
import pyproj
import shapely.geometry
from shapely.geometry import shape
import shapely.ops
import arcpy
import collections


def merge_mqual_geojson():
    shp_dir = r'X:\iocm_deliverables\ocs\GL07-16-TB-J\michigan\PRODUCT_CELLS'

    p_dir = Path(shp_dir)
    geojsons = p_dir.rglob('*_mqual.geojson')

    mqual_polys = None
    sr = arcpy.SpatialReference(26916)

    for i, g in enumerate(list(geojsons)):
        with open(g) as f:
            try:
                print(g)

                with open(g, 'r') as j:
                    mqual = json.load(j)['features']

                gc = [shape(poly["geometry"]).buffer(0) for poly in mqual]

                if i == 0:
                    mqual_polys = gc
                    properties = mqual[0]['properties']
                else:
                    mqual_polys = mqual_polys + gc
            except Exception as e:
                print(e)

    print('merging {} shapely polygons...'.format(len(mqual_polys)))
    merged_mqual = shapely.ops.cascaded_union(mqual_polys)

    print('copying features to shapefile...')
    mqual_shp_dir = r'X:\iocm_deliverables\ocs\GL07-16-TB-J\michigan\PRODUCT_CELLS'
    mqual_shp = r'X:\iocm_deliverables\ocs\GL07-16-TB-J\michigan\PRODUCT_CELLS\mqual.shp'

    fields = collections.OrderedDict([
        ('CATZOC', 3), 
        ('POSACC', 1), 
        ('SOUACC', 0.5), 
        ('SUREND', '20090327'),
        ('SURSTA', '20090327'), 
        ('TECSOU', 7), 
        ('INFORM', 'NOAA NGS-RSD'), 
        ('SORDAT', '20090327'), 
        ('SORIND', 'US,US,graph'), 
        ('FCSubtype', 40)
        ])

    arcpy.CreateFeatureclass_management(mqual_shp_dir, 'mqual.shp', 'POLYGON',
                                        str(g).replace('.geojson', '.shp'), 
                                        spatial_reference=sr)

    mquals = arcpy.da.InsertCursor(mqual_shp, ['SHAPE@WKT'] + list(fields.keys()))
    mquals.insertRow([merged_mqual.to_wkt()] + list(fields.values()))
    del mquals
