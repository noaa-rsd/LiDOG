""" 
LiDOG, the Lidar Deliverables for OCS Generator, processes 
RSD and/or JALBTCX bathymetric lidar DEMs to generate a generaized
product DEM and an associated M_QUAL shapefile.

Nick Forfinski-Sarkozi, NOAA Remote Sensing Division
nick.forfinski-sarkozi@noaa.gov
"""


import arcpy
import re
import xml.etree.ElementTree as ET
from arcpy.sa import *
import os
import collections
from pathlib import Path


class SourceDem:

    def __init__(self, path, lidog):
        self.path = path
        self.basename_unvalidated = self.path.stem
        self.basename = arcpy.ValidateTableName(self.basename_unvalidated)
        self.name = self.basename + '.tif'
        self.agg_factor = 5
        self.proj_dir = lidog.out_dir
        self.band4_cells = None
        self.raster = arcpy.sa.Raster(str(self.path))

    def aggregate(self, mosaic_path):
        src_res = '1m'
        prod_res = '5m' 
        arcpy.AddMessage('Aggregating mosaicked source DEM to product DEM...')
        agg_dem = Aggregate(str(mosaic_path), self.agg_factor, 'MINIMUM', 'TRUNCATE', 'NODATA')
        agg_dem_path = Path(str(mosaic_path).replace(src_res, prod_res))
        #agg_dem.save(str(agg_dem_path))
        return agg_dem, agg_dem_path

    def get_coverage_fc(self):
        arcpy.AddMessage('determing band-4 coverage of {}'.format(self.path))
        src_dem_coverage_fc = 'in_memory\src_dem_coverage'
        r = Aggregate(self.raster, 10, extent_handling='TRUNCATE')
        arcpy.RasterToPolygon_conversion(Int(r * 0), src_dem_coverage_fc, simplify='NO_SIMPLIFY')
        return src_dem_coverage_fc

    def extract_band4_cells(self):
        arcpy.MakeFeatureLayer_management(str(ProductCell.band4_shp), 'band4')
        extent_poly = self.get_coverage_fc()
        arcpy.SelectLayerByLocation_management('band4', 'INTERSECT', extent_poly)
        num_cells = int(arcpy.GetCount_management('band4').getOutput(0))
        arcpy.AddMessage('{} coverage intersects {} band-4 cells'.format(self.name, num_cells))
        self.band4_cells = self.proj_dir / (self.basename + '_BAND4_cells.shp')
        arcpy.CopyFeatures_management('band4', str(self.band4_cells))
        
        band4_cells = []
        with arcpy.da.SearchCursor(str(self.band4_cells), ['CellName', 'SHAPE@']) as cells:
            for cell in cells:
                band4_cells.append(cell)
        return band4_cells


class ProductDem:

    def __init__(self, agg_raster, agg_path):
        self.path = agg_path
        self.raster = agg_raster
        self.max_depth = -20

    def mask_land(self):
        query_str = 'VALUE >= 0 OR VALUE <= {}'.format(self.max_depth)
        agg_dem_water = SetNull(self.raster, self.raster, query_str)
        return agg_dem_water

    def generalize_water_coverage(self):
        arcpy.AddMessage('generalizing and land-masking product DEM...')
        generalized_dem5 = Aggregate(Int(self.mask_land()), 4, 'MINIMUM', 'TRUNCATE', 'NODATA')
        return generalized_dem5
        
    def trim(self, mqual_path):
        arcpy.AddMessage('clipping original product DEM with M_QUAL...')
        arcpy.Clip_management(self.raster, '#', str(self.path).replace('.tif', '_M_QUAL.tif'), 
                              str(mqual_path), clipping_geometry=True)


class ProductCell:

    band4_shp = Path('../support_shps/All_Band4_V5.shp')
    fields = ['SHAPE@', 'CellName']

    def __init__(self, cell, lidog, mqual):
        self.name = cell[0]
        self.geom = cell[1]
        self.proj_dir = lidog.out_dir
        self.product_cell_name = '{}_ProductCell'.format(self.name)
        self.mqual_name = None
        self.mqual_path = None
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
        arcpy.AddMessage('converting generalized DEM to polygons...')
        dem_polys = 'in_memory\polys'
        arcpy.RasterToPolygon_conversion(generalized_coverage, dem_polys, 'NO_SIMPLIFY')
        polys_dissolved = 'in_memory\DissolvedFC'
        arcpy.Dissolve_management(dem_polys, polys_dissolved, multi_part='SINGLE_PART')

        arcpy.AddMessage('denoising polygons...')
        self.mqual_name = '{}_M_QUAL.shp'.format(self.name)
        self.mqual_path = self.product_cell_path / self.mqual_name

        arcpy.CreateFeatureclass_management('in_memory', 'pre_mqual', 
                                            'POLYGON', str(self.mqual_schema_rpath), 
                                            spatial_reference=self.spatial_ref)

        fields = ['SHAPE@'] + list(self.mqual_fields.keys())
        mqual_values = list(self.mqual_fields.values())

        with arcpy.da.SearchCursor(polys_dissolved, ['SHAPE@', 'SHAPE@AREA']) as SearchCursor:
            insert_cursor = arcpy.da.InsertCursor('in_memory\pre_mqual', fields)
            for poly in SearchCursor:
                if poly[1] > self.denoise_threshold:
                    insert_cursor.insertRow([poly[0]] + mqual_values)
        
        arcpy.AddMessage('trimming pre_M_QUAL with cell geometry...')
        self.mqual_trimmed_path = Path(str(self.mqual_path) + '_TRIMMED.shp')
        arcpy.Clip_analysis('in_memory\pre_mqual', str(self.extent_fc_path), 
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

    def mosaic(self, clipped_src_dems):
        inputs = ';'.join(clipped_src_dems)
        arcpy.AddMessage('mosaicking clipped source DEMs...'.format(self.name))
        mosaic_name = '{}_1m_DEM_MOSAIC.tif'.format(self.name)
        arcpy.MosaicToNewRaster_management(inputs, str(self.product_cell_path), mosaic_name, 
                                           pixel_type='32_BIT_FLOAT', number_of_bands=1,
                                           mosaic_method='MINIMUM')
        return self.product_cell_path / mosaic_name

    def clip_src_dem(self, dem):
        arcpy.AddMessage('clipping {} with {} cell buffer...'.format(dem.name, self.name))
        clipped_dem_name = dem.name.replace('.tif', '_{}.tif'.format(self.name))
        clipped_dem_path = self.clipped_dem_dir / clipped_dem_name
        arcpy.Clip_management(str(dem.raster), '#', str(clipped_dem_path), 
                              str(self.buffer_fc_path), clipping_geometry=True)

    def create_product_src(self):
        clipped_src_dems = [str(c) for c in self.clipped_dem_dir.glob('*_{}.tif'.format(self.name))]
        if len(clipped_src_dems) > 1:
            product_src_path = self.mosaic(clipped_src_dems)
        elif len(clipped_src_dems) == 1:
            arcpy.AddMessage('only 1 clipped src DEM in cell; copying it to product DEM source...')
            product_src_name = '{}_5m_DEM_MOSAIC.tif'.format(self.name)
            product_src_path = self.product_cell_path / product_src_name
            arcpy.CopyRaster_management(clipped_src_dems[0], str(product_src_path))

        return product_src_path


class MetaData:

    def __init__(self, lidog):
        self.template_path = Path('../metadata_xml/LiDOG_metadata_template.xml')
        self.path = lidog.out_dir / (lidog.project_id + '_metadata.xml')
        self.srs_wkt = lidog.spatial_ref
        self.central_meridian = self.srs_wkt.centralMeridianInDegrees
        self.utm_zone = int((self.central_meridian + 180.0) / 6) + 1
        self.xml_root = None
        self.sursta = lidog.sursta
        self.surend = lidog.surend

        self.metadata = {
            'begdate': self.format_date(self.sursta),
            'enddate': self.format_date(self.surend),
            'proj_id': lidog.project_id,  # TODO: add to title tag 
            'longcm': self.central_meridian,
            'utmzone': self.utm_zone,
            }

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
        
    def update_metadata(self):
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
        self.mqual_schema_rpath = Path('../support_shps/M_QUAL.shp')
        self.mquals = []
        self.cells = None
        self.proj_dir = lidog.out_dir

    def combine_mquals(self):
        arcpy.AddMessage('merging cell M-QUALs to create project-wide M_QUAL...')
        if len(self.mquals) > 1:
            temp_mqual = r'in_memory\temp_mqual_shp'
            current_mqual = r'in_memory\current_mqual'
            arcpy.Union_analysis([str(self.mquals[0]), str(self.mquals[1])], current_mqual)

            for i in range(2, len(self.mquals)):
                arcpy.Union_analysis([current_mqual, str(self.mquals[i])], temp_mqual)
                arcpy.CopyFeatures_management(temp_mqual, current_mqual)
            arcpy.Dissolve_management(current_mqual, str(self.project_mqual_path).replace('.shp', ''))

        else:
            arcpy.FeatureClassToFeatureClass_conversion(str(self.mquals[0]), str(self.proj_dir), 
                                                        self.project_mqual_name)


class LiDOG:
       
    def __init__(self):
        self.project_id = arcpy.GetParameterAsText(0)
        self.sursta = arcpy.GetParameterAsText(1)
        self.surend = arcpy.GetParameterAsText(2)
        self.spatial_ref = arcpy.GetParameter(3)
        self.source_dems = [Path(str(dem1)) for dem1 in arcpy.GetParameter(4)]
        self.z_convention = arcpy.GetParameterAsText(5)
        self.out_dir = Path(arcpy.GetParameterAsText(6))
        self.num_dems = len(self.source_dems)
        self.product_cells = {}


if __name__ == '__main__':

    lidog = LiDOG()
    
    metadata = MetaData(lidog)

    mqual = Mqual(lidog, metadata)

    arcpy.CheckOutExtension('Spatial')
    arcpy.env.workspace = str(lidog.out_dir)
    #arcpy.env.overwriteOutput = True

    for i, dem_path in enumerate(lidog.source_dems, 1):
        arcpy.AddMessage('{} (source DEM {} of {})'.format('*' * 40, i, lidog.num_dems))

        src_dem = SourceDem(dem_path, lidog)
        cells = src_dem.extract_band4_cells()

        num_dem_cells = len(cells)
        for j, cell in enumerate(cells, 1):
            cell_name = cell[0]
            arcpy.AddMessage('{} (DEM cell {} of {})'.format('-' * 40, j, num_dem_cells))
            arcpy.AddMessage('processing band-4 product cell {}...'.format(cell_name))

            if cell_name not in lidog.product_cells.keys():
                product_cell = ProductCell(cell, lidog, mqual)
                lidog.product_cells.update({cell_name: product_cell})

            lidog.product_cells[cell_name].clip_src_dem(src_dem)
        
    num_cells = len(lidog.product_cells)
    for i, (cell_name, cell) in enumerate(lidog.product_cells.items(), 1):     
        arcpy.AddMessage('{} cell {} ({} of {})'.format('=' * 40, cell_name, i, num_cells))

        product_src_path = cell.create_product_src()
        agg_dem, agg_dem_path = src_dem.aggregate(product_src_path)
        prod_dem = ProductDem(agg_dem, agg_dem_path)

        cell_mqual = cell.create_mqual(prod_dem.generalize_water_coverage())
        mqual.mquals.append(cell_mqual)
        prod_dem.trim(cell_mqual)

    arcpy.AddMessage('+' * 40)
    mqual.combine_mquals()

    arcpy.AddMessage('exporting xml metatdata file...')
    metadata.get_xml_template()
    metadata.update_metadata()
    metadata.export_metadata()
