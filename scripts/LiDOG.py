"""
This tool processes NOAA RSD and/or USACE JALBTCX DEM data 
and creates an M_QUAL polygon with appropriate attribution.

Originally developed based on "LidarProcessor" tool created 
by Noel Dyer.

Author:
Nick Forfinski-Sarkozi, NOAA Remote Sensing Division
nick.forfinski-sarkozi@noaa.gov
"""

import arcpy
import re
import xml.etree.ElementTree as ET
from arcpy.sa import *
import os
import json
import collections
from pathlib import Path
import rasterio
from rasterio.mask import mask
from shapely.geometry import shape, GeometryCollection
from shapely.ops import transform
from functools import partial
import pyproj
from rasterio.io import MemoryFile 
from rasterio import features
from rasterio.enums import Resampling
from rasterio import Affine
import numpy as np
import time


class SourceDem:

    def __init__(self, path, lidog):
        self.path = path
        self.basename_unvalidated = self.path.stem
        self.basename = arcpy.ValidateTableName(self.basename_unvalidated)
        self.name = self.basename + '.tif'
        self.agg_factor = 5
        self.proj_dir = lidog.out_dir
        self.band4_cells = None
        #self.raster = arcpy.sa.Raster(str(self.path))

    def aggregate(self, mosaic_path):
        arcpy.AddMessage('Aggregating mosaicked source DEM to product-DEM resolution...')
        agg_dem = Aggregate(str(mosaic_path), self.agg_factor, 'MINIMUM', 'TRUNCATE', 'NODATA')
        return agg_dem

    def get_coverage_fc(self):
        arcpy.AddMessage('determining coverage of DEM...')
        with rasterio.open(self.path) as r:
            t=r.meta['transform']
            data = r.read_masks(1, out_shape=(r.height // 10, r.width // 10), 
                                resampling=Resampling.average)

        # adjust transform after raster resampling
        transform = Affine(t.a * 10, t.b, t.c, t.d, t.e * 10, t.f)
        data = np.ma.array(data, mask=(data == 0))
        mask = features.shapes(data, mask=None, transform=transform)

        coverage_shp = 'in_memory\cov'
        sr = arcpy.SpatialReference(r.crs.to_epsg())
        arcpy.CreateFeatureclass_management('in_memory', 'cov', spatial_reference=sr)
        cursor = arcpy.da.InsertCursor(coverage_shp, ['SHAPE@WKT'])

        for poly, value in mask:
            poly_wkt = shape(poly).wkt
            cursor.insertRow([poly_wkt])
        del cursor
        return coverage_shp

    def extract_band4_cells(self):

        self.band4_cells_dir = self.proj_dir / 'Source_DEM_Band4_Cells'
        try:
            os.mkdir(self.band4_cells_dir)
        except Exception as e:
            pass

        arcpy.MakeFeatureLayer_management(str(ProductCell.band4_shp), 'band4')
        extent_poly = self.get_coverage_fc()
        arcpy.SelectLayerByLocation_management('band4', 'INTERSECT', extent_poly)
        num_cells = int(arcpy.GetCount_management('band4').getOutput(0))
        arcpy.AddMessage('{} coverage intersects {} band-4 cells'.format(self.name, num_cells))

        self.band4_cells = self.band4_cells_dir / (self.basename + '_BAND4_cells.shp')
        arcpy.CopyFeatures_management('band4', str(self.band4_cells))
        
        band4_cells = []
        with arcpy.da.SearchCursor(str(self.band4_cells), ['CellName', 'SHAPE@']) as cells:
            for cell in cells:
                band4_cells.append(cell)
        return band4_cells, self.band4_cells


class ProductDem:

    def __init__(self, agg_raster, cell):
        self.product_cell_name = cell.name
        self.product_cell_path = cell.product_cell_path
        self.pre_product_path = self.product_cell_path / (self.product_cell_name + '_5m_DEM.tif')
        self.raster = agg_raster
        self.raster.save(str(self.pre_product_path))
        #self.max_depth = -20

    def mask_land(self):
        #query_str = 'VALUE >= 0 OR VALUE <= {}'.format(self.max_depth)
        query_str = 'VALUE >= 0'
        agg_dem_water = SetNull(self.raster, self.raster, query_str)
        return agg_dem_water

    def generalize_water_coverage(self):
        arcpy.AddMessage('generalizing preliminary product DEM water coverage...')
        generalized_dem5 = Aggregate(Int(self.mask_land()), 4, 'MINIMUM', 'TRUNCATE', 'NODATA')
        return generalized_dem5
        

class ProductCell:

    band4_shp = Path('../support_files/All_Band4_V5.shp')
    src_dems_template_shp = Path(r'../support_files/src_dem_extent_shp_template.shp')
    fields = ['SHAPE@', 'CellName']

    def __init__(self, cell, lidog, mqual):
        self.name = cell[0]
        self.geom = cell[1]
        self.proj_dir = lidog.out_dir
        self.product_cell_name = '{}_ProductCell'.format(self.name)
        self.mqual_name = None
        self.mqual_path = None
        self.mqual_in_memory = None
        self.mqual_trimmed_path = None
        self.mqual_fields = mqual.fields
        self.mqual_schema_rpath = mqual.mqual_schema_rpath
        self.spatial_ref = lidog.spatial_ref
        self.buffer_meters = 50
        self.denoise_threshold = 1200  # square meters
        self.product_cell_path = self.make_cell_dir()
        self.clipped_dem_dir = self.product_cell_path / 'source_dems'
        self.extent_fc_path = self.make_extent_fc()
        self.buffer_fc_path = self.buffer()

        try:
            os.mkdir(str(self.clipped_dem_dir))
        except Exception as e:
            pass

        pass

    def create_mqual(self, generalized_coverage):
        arcpy.AddMessage('generating M_QUAL from generalized water coverage...')
        dem_polys = 'in_memory\polys'
        arcpy.RasterToPolygon_conversion(generalized_coverage, dem_polys, 'NO_SIMPLIFY')
        polys_dissolved = 'in_memory\DissolvedFC'
        arcpy.Dissolve_management(dem_polys, polys_dissolved, multi_part='SINGLE_PART')

        self.mqual_name = '{}_M_QUAL.shp'.format(self.name)
        self.mqual_path = self.product_cell_path / self.mqual_name

        self.mqual_in_memory ='in_memory\pre_mqual_{}'.format(cell.name)
        arcpy.CreateFeatureclass_management('in_memory', 'pre_mqual_{}'.format(cell.name), 
                                            'POLYGON', str(self.mqual_schema_rpath), 
                                            spatial_reference=self.spatial_ref)

        fields = ['SHAPE@'] + list(self.mqual_fields.keys())
        mqual_values = list(self.mqual_fields.values())

        with arcpy.da.SearchCursor(polys_dissolved, ['SHAPE@', 'SHAPE@AREA']) as SearchCursor:
            insert_cursor = arcpy.da.InsertCursor(self.mqual_in_memory, fields)
            for poly in SearchCursor:
                if poly[1] > self.denoise_threshold:
                    insert_cursor.insertRow([poly[0]] + mqual_values)
        
        self.mqual_trimmed_path = Path(str(self.mqual_path) + '_TRIMMED.shp')
        arcpy.Clip_analysis(self.mqual_in_memory, str(self.extent_fc_path), 
                            str(self.mqual_path))
                    
        return self.mqual_path

    def make_cell_dir(self):
        cell_path = self.proj_dir / self.product_cell_name
        if not arcpy.Exists(str(cell_path)):
            try:
                os.mkdir(str(cell_path))
                return cell_path
            except Exception as e:
                arcpy.AddMessage(e)
        else:
            return cell_path

    def make_extent_fc(self):
        extent_fc_path = self.product_cell_path / (self.name + '_CellOutline.shp')
        if not arcpy.Exists(str(extent_fc_path)):
            try:
                arcpy.CopyFeatures_management(self.geom, str(extent_fc_path))
                return extent_fc_path
            except Exception as e:
                arcpy.AddMessage(e)
        else:
            return extent_fc_path

    def buffer(self):
        buffer_fc_path = r'in_memory\buffer_{}'.format(self.name)
        if not arcpy.Exists(str(buffer_fc_path)):
            try:
                arcpy.Buffer_analysis(str(self.extent_fc_path), 
                                      str(buffer_fc_path), 
                                      '{} Meters'.format(self.buffer_meters))
                return buffer_fc_path
            except Exception as e:
                arcpy.AddMessage(e)            
        else:
            return buffer_fc_path

    def mosaic(self):
        clipped_src_dems = [str(c) for c in self.clipped_dem_dir.glob('*_{}.tif'.format(self.name))]
        inputs = ';'.join(clipped_src_dems)
        arcpy.AddMessage('mosaicking clipped source DEMs...'.format(self.name))
        mosaic_name = '{}_src'.format(self.name)
        arcpy.MosaicToNewRaster_management(inputs, 'in_memory', mosaic_name, 
                                           pixel_type='32_BIT_FLOAT', number_of_bands=1, 
                                           mosaic_method='MINIMUM')

        mosaic_path = 'in_memory\{}'.format(mosaic_name)
        mosaic_resolution = arcpy.Describe(mosaic_path).meancellwidth

        return mosaic_path, mosaic_resolution

    def get_cell_geometry(self):
        geojson_path = str(self.extent_fc_path).replace('.shp', '.geojson')
        arcpy.FeaturesToJSON_conversion(self.buffer_fc_path, geojson_path, geoJSON='GEOJSON')
        with open(geojson_path, 'r') as j:
            cell_poly = json.load(j)['features'][0]['geometry']
        cell_geom = shape(cell_poly)  # convert into shapely geometry
        return cell_geom

    def get_mqual_geometry(self):
        geojson_path = str(self.mqual_path).replace('.shp', '.geojson')
        arcpy.FeaturesToJSON_conversion(self.mqual_in_memory, geojson_path, geoJSON='GEOJSON')
        with open(geojson_path, 'r') as j:
            mqual = json.load(j)['features']
        gc = GeometryCollection([shape(poly["geometry"]) for poly in mqual])
        return gc

    def mask_dem(self, dem_path, geom, masked_dem_path):
        src_r = rasterio.open(dem_path)
        project = partial(
            pyproj.transform,
            pyproj.Proj(init='epsg:4326'),
            pyproj.Proj(init='epsg:{}'.format(src_r.crs.to_epsg())))

        geom_utm = transform(project, geom)
        if isinstance(geom_utm, GeometryCollection):
            geom = geom_utm
        else:  # individual shapely polygons
            geom = [geom_utm]

        try:
            out_r, out_transform = mask(dataset=src_r, shapes=geom, crop=True)
            out_meta = src_r.meta.copy()
            out_meta.update({'driver': 'GTiff',
                             'height': out_r.shape[1],
                             'width': out_r.shape[2],
                             'transform': out_transform,
                             'crs': src_r.crs.to_proj4(),
                             'compress': 'lzw'})

            with rasterio.open(masked_dem_path, "w", **out_meta) as masked_dem:
                masked_dem.write(out_r)

        except Exception as e:
            arcpy.AddMessage(e)

    def clip_src_dem(self, dem):
        arcpy.AddMessage('clipping {} with {} cell buffer...'.format(dem.name, self.name))
        clipped_dem_name = dem.name.replace('.tif', '_{}.tif'.format(self.name))
        clipped_dem_path = self.clipped_dem_dir / clipped_dem_name
        cell_geom = self.get_cell_geometry()
        self.mask_dem(dem.path, cell_geom, clipped_dem_path)

    def clip_pre_product_dem(self, dem):
        arcpy.AddMessage('clipping preliminary product DEM with M_QUAL...')
        product_dem_path = self.product_cell_path / (self.product_cell_name + '_5m_DEM_M_QUAL.tif')
        mqual_geom = self.get_mqual_geometry()
        self.mask_dem(dem, mqual_geom, product_dem_path)


class MetaData:

    def __init__(self, lidog):
        self.template_path = Path('../metadata_xml/LiDOG_metadata_template.xml')
        self.meta_library_path = Path('../support_files/meta_library.json')
        self.path = lidog.out_dir / (lidog.project_id + '_metadata.xml')
        self.srs_wkt = lidog.spatial_ref
        self.central_meridian = self.srs_wkt.centralMeridianInDegrees
        self.utm_zone = int((self.central_meridian + 180.0) / 6) + 1
        self.xml_root = None
        self.sursta = lidog.sursta
        self.surend = lidog.surend
        self.data_src = lidog.data_src
        self.meta_library = self.get_meta_library()
        self.bounding_coordinates = None

        self.metadata = {
            'begdate': self.format_date(self.sursta),
            'enddate': self.format_date(self.surend),
            'proj_id': lidog.project_id,  # TODO: add to title tag 
            'longcm': self.central_meridian,
            'utmzone': self.utm_zone,
            'procdesc': self.meta_library['procdesc'][self.data_src],
            'westbc': None,
            'eastbc': None,
            'northbc': None,
            'southbc': None,
            }

    def popuate_extents(self, mqual_path):
        mqual_extents = arcpy.Describe(str(mqual_path)).extent
        mqual_extents_DD = mqual_extents.projectAs(arcpy.SpatialReference(4326))
        self.metadata['westbc'] = mqual_extents_DD.XMin
        self.metadata['eastbc'] = mqual_extents_DD.XMax
        self.metadata['northbc'] = mqual_extents_DD.YMin
        self.metadata['southbc'] = mqual_extents_DD.YMax

    def get_meta_library(self):
        lidog_dir = os.path.dirname(os.path.realpath(__file__))
        os.chdir(lidog_dir)
        with open(self.meta_library_path) as json_file:
            meta_library = json.load(json_file)
        return meta_library

    @staticmethod
    def format_date(date):
        """reformats date from mm/dd/yyyy to yyyymmdd"""

        d = date.split('/')
        return d[-1] + d[0] + d[1]

    def get_xml_template(self):
        lidog_dir = os.path.dirname(os.path.realpath(__file__))
        os.chdir(lidog_dir)
        tree = ET.parse(self.template_path)
        self.xml_root = tree.getroot()
        
    def update_metadata(self, mqual_path):
        self.popuate_extents(mqual_path)
        for metadatum, val in self.metadata.items():
            for e in self.xml_root.iter(metadatum):
                e.text = str(val)
    
    def export_metadata(self):
        updated_xml = ET.tostring(self.xml_root).decode('utf-8')
        with open(self.path, "w") as f:
            f.write(updated_xml)
        
        pass


class Mqual:

    def __init__(self, lidog, meta):
        self.project_mqual_name = '{}_M_QUAL.shp'.format(lidog.project_id)
        self.project_mqual_path = lidog.out_dir / self.project_mqual_name
        self.dem_mqual_path = None
        self.fields = collections.OrderedDict([
            ('CATZOC', 3), 
            ('POSACC', 1), 
            ('SOUACC', 0.5), 
            ('SUREND', meta.surend),
            ('SURSTA', meta.sursta), 
            ('TECSOU', 7), 
            ('INFORM', 'NOAA NGS-RSD'), 
            ('SORDAT', lidog.surend), 
            ('SORIND', 'US,US,graph'), 
            ('FCSubtype', 40)
            ])
        self.sursta = lidog.sursta
        self.surend = lidog.surend
        self.spatial_ref = lidog.spatial_ref
        self.mqual_schema_rpath = Path('../support_files/M_QUAL_TEMPLATE.shp')
        self.mquals = []
        self.cells = None
        self.proj_dir = lidog.out_dir

    def combine_mquals(self):
        arcpy.AddMessage('merging cell M-QUALs to create project-wide M_QUAL...')
        if len(self.mquals) > 1:
            temp_mqual = r'in_memory\temp_mqual_shp'
            current_mqual = r'in_memory\current_mqual'
            arcpy.Union_analysis([str(self.mquals[0]), 
                                  str(self.mquals[1])], 
                                 current_mqual)

            for i in range(2, len(self.mquals)):
                arcpy.Union_analysis([current_mqual, str(self.mquals[i])], temp_mqual)
                arcpy.CopyFeatures_management(temp_mqual, current_mqual)
            arcpy.Dissolve_management(current_mqual, str(self.project_mqual_path))

        else:
            arcpy.FeatureClassToFeatureClass_conversion(str(self.mquals[0]), 
                                                        str(self.proj_dir), 
                                                        self.project_mqual_name)

        return self.project_mqual_path


class LiDOG:
       
    def __init__(self):
        self.project_id = arcpy.GetParameterAsText(0)
        self.data_src = arcpy.GetParameterAsText(1)
        self.sursta = arcpy.GetParameterAsText(2)
        self.surend = arcpy.GetParameterAsText(3)
        self.spatial_ref = arcpy.GetParameter(4)
        self.source_dems = [Path(str(dem1)) for dem1 in arcpy.GetParameter(5)]
        self.z_convention = arcpy.GetParameterAsText(6)
        self.out_dir = Path(arcpy.GetParameterAsText(7))
        self.num_dems = len(self.source_dems)
        self.product_cells = {}
        self.src_dem_band4_cells = []
        self.src_dem_band4_shp_dir = self.out_dir
        self.project_band4_cells_path = None
        self.src_dems_extent_path = None
        self.mqual_path = None

        try:
            os.mkdir(str(self.src_dem_band4_shp_dir))
        except Exception as e:
            pass

    def create_project_band4_cells_shp(self):
        project_band4_cells_name = self.project_id + '_Band4_Cells.shp'
        self.project_band4_cells_path = self.src_dem_band4_shp_dir / project_band4_cells_name
        proj_band4_shp_temp = r'in_memory\band4'
        arcpy.CreateFeatureclass_management('in_memory', 'band4', 'POLYGON',
                                            str(ProductCell.band4_shp),
                                            spatial_reference=self.spatial_ref)
        band4_cells = arcpy.da.InsertCursor(proj_band4_shp_temp, ['CellName', 'SHAPE@'])
        for cell in self.product_cells.values(): 
            band4_cells.insertRow([cell.name, cell.geom])
        del band4_cells
        arcpy.CopyFeatures_management(proj_band4_shp_temp, str(self.project_band4_cells_path))

    def create_source_dems_extents(self):
        src_dems_extent_name = self.project_id + '_Source_Dem_Extents.shp'
        self.src_dems_extent_path = self.out_dir / src_dems_extent_name
        src_dems_extent_name_temp = r'in_memory\extents'
        arcpy.CreateFeatureclass_management('in_memory', 'extents', 'POLYGON',
                                            str(ProductCell.src_dems_template_shp),
                                            spatial_reference=self.spatial_ref)
        src_dem_extents = arcpy.da.InsertCursor(src_dems_extent_name_temp, ['name', 'resolution', 'SHAPE@'])
        for dem in self.source_dems: 
            desc = arcpy.Describe(str(dem))
            dem_res = desc.meancellwidth
            src_dem_extents.insertRow([dem.name, dem_res, desc.extent.polygon])
        del src_dem_extents
        arcpy.CopyFeatures_management(src_dems_extent_name_temp, str(self.src_dems_extent_path))       

    def generate_summary_plot(self):
        arcpy.AddMessage('creating ArcGIS project with results...')
        project_aprx = self.out_dir / (self.project_id + '.aprx')

        try:
            os.mkdir(str(aprx_dir))
        except Exception as e:
            pass

        def buffer_extent(ext, buffer_factor):
            dx = ext.XMax - ext.XMin
            dy = ext.YMax - ext.YMin
            buff_x = dx * buffer_factor
            buff_y = dy * buffer_factor 
            ext.XMin -= buff_x
            ext.XMax += buff_x
            ext.YMin -= buff_y
            ext.YMax += buff_y
            return ext, dx, dy
        
        def determine_layout_orientation(dx, dy):
            if dy >= dx:
                orientation = 'portrait'
            else:
                orientation = 'landscape'
            return orientation

        buffer_factor = 0.1
        ext = arcpy.Describe(str(self.project_band4_cells_path)).extent
        ext, dx, dy = buffer_extent(ext, buffer_factor)
        orientation = determine_layout_orientation(dx, dy)

        lidog_dir = Path(os.path.dirname(os.path.realpath(__file__)))
        support_files_dir = lidog_dir.parent / 'support_files'

        if orientation == 'portrait':
            aprx_template_path = lidog_dir.parent / 'LiDOG_ProjectWideTemplate_PORTRAIT.aprx'
        elif orientation == 'landscape':
            aprx_template_path = lidog_dir.parent / 'LiDOG_ProjectWideTemplate_LANDSCAPE.aprx'

        lyr_cells_template_path = support_files_dir / 'Project_Band4_Cells_TEMPLATE.lyrx'
        lyr_extents_template_path = support_files_dir / 'Source_Dem_Extents_TEMPLATE.lyrx'
        lyr_mqual_template_path = support_files_dir / 'M_QUAL_TEMPLATE.lyrx'

        aprx = arcpy.mp.ArcGISProject(str(aprx_template_path))
        aprx_map = aprx.listMaps("Map")[0]
        
        lyrx_cells = self.project_band4_cells_path.parent / (self.project_band4_cells_path.stem + '.lyrx')
        lyrx_extents = self.src_dems_extent_path.parent / (self.src_dems_extent_path.stem + '.lyrx')
        lyrx_mqual = self.mqual_path.parent / (self.mqual_path.stem + '.lyrx')

        arcpy.SaveToLayerFile_management(arcpy.MakeFeatureLayer_management(str(self.project_band4_cells_path)), str(lyrx_cells))
        arcpy.SaveToLayerFile_management(arcpy.MakeFeatureLayer_management(str(self.src_dems_extent_path)), str(lyrx_extents))
        arcpy.SaveToLayerFile_management(arcpy.MakeFeatureLayer_management(str(self.mqual_path)), str(lyrx_mqual))

        lyr_cells_template = arcpy.mp.LayerFile(str(lyr_cells_template_path))
        lyr_extents_template = arcpy.mp.LayerFile(str(lyr_extents_template_path))
        lyr_mqual_template = arcpy.mp.LayerFile(str(lyr_mqual_template_path))

        lyrx_cells = arcpy.mp.LayerFile(str(lyrx_cells))
        lyrx_extents = arcpy.mp.LayerFile(str(lyrx_extents))
        lyrx_mqual = arcpy.mp.LayerFile(str(lyrx_mqual))

        arcpy.ApplySymbologyFromLayer_management(lyrx_cells.listLayers()[0], lyr_cells_template.listLayers()[0])
        arcpy.ApplySymbologyFromLayer_management(lyrx_extents.listLayers()[0], lyr_extents_template.listLayers()[0])
        arcpy.ApplySymbologyFromLayer_management(lyrx_mqual.listLayers()[0], lyr_mqual_template.listLayers()[0])

        lyrx_cells.save()
        lyrx_extents.save()
        lyrx_mqual.save()

        aprx_map.addLayer(lyrx_cells, 'BOTTOM')
        aprx_map.addLayer(lyrx_mqual, 'BOTTOM')
        aprx_map.addLayer(lyrx_extents, 'BOTTOM')

        lyrx_cells_lyr_basename = lyrx_cells.listLayers()[0].name
        lyr = aprx_map.listLayers(lyrx_cells_lyr_basename)[0]
        lyt = aprx.listLayouts("Layout")[0]
        mf = lyt.listElements("mapframe_element", "Map Frame")[0]

        lyt.listElements("text_element", "Text")[0].text = 'Band-4 Product Cell Index Map\n{}'.format(self.project_id)

        ext = mf.getLayerExtent(lyr, False, True)
        ext, __, __ = buffer_extent(ext, buffer_factor)

        mf.camera.setExtent(ext)

        aprx.saveACopy(str(project_aprx))

        arcpy.AddMessage('generating summary plot pdf...')
        pdf_path = str(self.out_dir / '{}_Product_Cells_Overview.pdf'.format(self.project_id))

        arcpy.AddMessage('generating summary plot png...')
        png_path = str(self.out_dir / '{}_Product_Cells_Overview.png'.format(self.project_id))

        lyt.exportToPDF(pdf_path, resolution=600)
        lyt.exportToPNG(png_path, resolution=600)


def set_env_vars():
    user_dir = os.path.expanduser('~')

    script_path = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                          'anaconda3', 'Scripts')

    gdal_data = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                        'anaconda3', 'envs', 'LiDOG', 
                                        'Library', 'share', 'gdal')

    proj_lib = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                       'anaconda3', 'envs', 'LiDOG', 
                                       'Library', 'share')

    if script_path.name not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep + str(script_path)

    os.environ["GDAL_DATA"] = str(gdal_data)
    os.environ["PROJ_LIB"] = str(proj_lib)


if __name__ == '__main__':
    set_env_vars()

    lidog = LiDOG()
    metadata = MetaData(lidog)
    mqual = Mqual(lidog, metadata)

    arcpy.CheckOutExtension('Spatial')
    arcpy.env.workspace = str(lidog.out_dir)
    arcpy.env.overwriteOutput = True

    for i, dem_path in enumerate(lidog.source_dems, 1):
        arcpy.AddMessage('{} (source DEM {} of {})'.format('*' * 60, i, lidog.num_dems))

        # clip 1-m DEM to band 4 cells
        src_dem = SourceDem(dem_path, lidog)
        cells, cells_shp = src_dem.extract_band4_cells()
        #lidog.src_dem_band4_cells.append(cells_shp)

        # loop through each cell intersecting dem1
        num_dem_cells = len(cells)
        for j, cell in enumerate(cells, 1):
            cell_name = cell[0]
            arcpy.AddMessage('{} (DEM cell {} of {})'.format('-' * 40, j, num_dem_cells))
            arcpy.AddMessage('processing band-4 product cell {}...'.format(cell_name))

            # create product cell (if not already created)
            if cell_name not in lidog.product_cells.keys():
                product_cell = ProductCell(cell, lidog, mqual)
                lidog.product_cells.update({cell_name: product_cell})

            # clip dem1 with cell buffer
            lidog.product_cells[cell_name].clip_src_dem(src_dem)
        
    # create shapefile with source DEM extents
    lidog.create_source_dems_extents()

    num_cells = len(lidog.product_cells)
    for i, (cell_name, cell) in enumerate(lidog.product_cells.items(), 1):     
        arcpy.AddMessage('{} cell {} ({} of {})'.format('=' * 60, cell_name, i, num_cells))
        src_dem_mosaic_path, src_res = cell.mosaic()

        if src_res == 1:
            product_source_dem = src_dem.aggregate(src_dem_mosaic_path)
        else:
            arcpy.AddMessage('source DEM resolution != 1m, not aggregating source DEM')
            product_source_dem = arcpy.Raster(str(src_dem_mosaic_path))
        
        prod_dem = ProductDem(product_source_dem, cell)

        # create cell mqual
        cell_mqual = cell.create_mqual(prod_dem.generalize_water_coverage())
        mqual.mquals.append(cell_mqual)

        # clip pre product DEM with mqual
        cell.clip_pre_product_dem(prod_dem.pre_product_path)

    # create project-wide M_QUAL
    arcpy.AddMessage('+' * 60)
    lidog.mqual_path = mqual.combine_mquals()

    # create project-wide band4 cells shp
    arcpy.AddMessage('creating project-wide band4 cells shp...')
    lidog.create_project_band4_cells_shp()

    # export pdf of product cells overview
    lidog.generate_summary_plot()
    
    # update metadata with project-wide M_QUAL extents
    arcpy.AddMessage('exporting xml metatdata file...')
    metadata.get_xml_template()
    metadata.update_metadata(lidog.mqual_path)
    metadata.export_metadata()
