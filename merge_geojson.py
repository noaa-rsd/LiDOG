from pathlib import Path
from functools import partial
import json
import pyproj
import shapely.geometry
from shapely.geometry import shape
import shapely.ops
import arcpy


shp_dir = r'X:\iocm_deliverables\ocs\GL07-16-TB-J\michigan\PRODUCT_CELLS'
from pathlib import Path
p_dir = Path(shp_dir)
geojsons = p_dir.rglob('*_mqual.geojson')

mqual_polys = None

for i, g in enumerate(list(geojsons)):
    with open(g) as f:
        try:
            print(g)

            with open(g, 'r') as j:
                mqual = json.load(j)['features']
                print(json.dumps(mqual[0], indent=1))
            gc = [shape(poly["geometry"]).buffer(0) for poly in mqual]

            if i == 0:
                mqual_polys = gc
            else:
                mqual_polys = mqual_polys + gc
        except Exception as e:
            print(e)


#print('merging {} shapely polygons...'.format(len(mqual_polys)))
#merged_mqual = shapely.ops.cascaded_union(mqual_polys)
#print(merged_mqual)
#arc_polys = arcpy.FromWKT(merged_mqual.to_wkt(), arcpy.SpatialReference(26916))
#print(arc_polys)

#mqual_shp = r'X:\iocm_deliverables\ocs\GL07-16-TB-J\michigan\PRODUCT_CELLS\mqual.shp'
#arcpy.CopyFeatures_management(arc_polys, mqual_shp)



