import os
import json
import math
from pathlib import Path
import subprocess
import numpy as np
import rasterio
from rasterio import features
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
import pathos.pools as pp


class VerticalTransform:

    UTM = {'init': 'epsg:6346'}

    def __init__(self, dirs, params):
        self.dirs = dirs
        self.params = params
        self.vdatum_gdf = gpd.read_file(self.dirs['mcu_regions_path']).to_crs(self.UTM)
        self.xyz = None
        self.chunk_size = 1000000
        self.num_chunks = None

    def mask_dem(self, dem, geom, masked_dem_path):
        out_crs = 'epsg:{}'.format(self.spatial_ref.PCSCode)

        # only tranform cell boundary, not mqual
        if  not isinstance(geom, GeometryCollection):
            project = partial(
                pyproj.transform,
                pyproj.Proj(init='epsg:4326'),
                pyproj.Proj(init=out_crs))
            geom = [transform(project, geom)]

        try:
            dem_masked, transform = mask(dataset=dem, shapes=geom, crop=True)
            meta = dem.meta.copy()
            meta.update({
                'driver': 'GTiff',
                'height': dem_masked.shape[1],
                'width': dem_masked.shape[2],
                'transform': transform,
                'count': 1,
                'crs': out_crs,
                'compress': 'lzw'
                })

            with rasterio.open(masked_dem_path, "w", **meta) as masked_dem:
                masked_dem.write(dem_masked[0], 1)

        except Exception as e:
            arcpy.AddMessage(e)

    def merge_sep_models(self, dem):
        nodata_val = -999999
        sep_paths = list(self.dirs['vdatum_out_dir'].glob('*.txt'))
        df = pd.concat((pd.read_csv(f, skiprows=0) for f in sep_paths))
        df = df.sort_values(by=['1', '0'])
        null_model = df.to_numpy()[:, 2].reshape(dem.shape)
        null_model = np.where(null_model == nodata_val, np.nan, null_model)
        print(null_model)
        with rasterio.open(
                self.dirs['vdatum_out_dir'] / 'MarcoIsland_SEPARATION.tif', 'w',
                driver='GTiff', 
                dtype=rasterio.float64, 
                count=1, 
                width=dem.width, 
                height=dem.height, 
                transform=dem.transform,
                nodata=np.nan) as dst:
            dst.write(np.flipud(null_model), indexes=1)

    def create_null_tiff(self, dem, null_path):
        print('creating null tiff...')
        profile = dem.profile
        mask = dem.read() == profile['nodata']
        mask = mask.astype('ubyte').squeeze()
        print(mask)
        print(mask.dtype)
        profile.update({'count': 1, 
                        'dtype': rasterio.ubyte,
                        'nodata': 1})
        with rasterio.open(null_path, 'w', **profile) as dst:
            dst.write(mask, indexes=1)

    def create_null_ascii(self, dem):
        print('creating null ascii...')
        elevation = dem.read().squeeze()
        t = dem.profile['transform']
        eastings = t.c
        northings = t.f
        x = np.arange(0, dem.shape[1]) * t.a + t.c 
        y = np.arange(0, dem.shape[0]) * t.e + t.f
        xs, ys = np.meshgrid(x, y)
        xy = np.vstack((xs.ravel(), ys.ravel())).T
        self.xyz = np.hstack((xy, np.zeros((xy.shape[0], 1))))
        self.num_chunks = math.ceil(self.xyz.shape[0] / self.chunk_size)

    def create_sep_model(self, null_path):
        import os
        from pathlib import Path
        import subprocess

        print('creating separation model...')
        ascii_files = f'{null_path};{self.dirs["vdatum_out_dir"]}'
        os.chdir(self.dirs['jar_dir'])

        command_str = f'{str(self.dirs["java_dir"])} -jar vdatum.jar '
        command_str += f'ihorz:{self.params["ihorz"]} '
        command_str += f'ivert:{self.params["ivert"]} '
        command_str += f'ohorz:{self.params["ohorz"]} '
        command_str += f'overt:{self.params["overt"]} '
        command_str += f'-file:txt:comma,0,1,2,skip1:{ascii_files} '
        command_str += 'region:3'
        print(command_str)

        orig_PATH = Path(os.environ["PATH"])
        orig_GDAL_DATA = Path(os.environ["GDAL_DATA"])

        PATH = str(self.dirs['gdal_dir'])
        PATH += os.pathsep + str(self.dirs['gdal_dir'] / 'bin')
        PATH += os.pathsep + str(self.dirs['gdal_dir'].joinpath('bin', 'gdal', 'aps'))
        GDAL_DATA = str(self.dirs['gdal_dir'].joinpath('bin', 'gdal-data'))
        GDAL_DRIVER_PATH = str(self.dirs['gdal_dir'].joinpath('bin', 'gdal', 'plugins'))
        os.environ["PATH"] = PATH
        os.environ["GDAL_DATA"] = GDAL_DATA
        os.environ["GDAL_DRIVER_PATH"] = GDAL_DRIVER_PATH

        try:
            print('running VDatum...')
            subprocess.call(command_str.split(' '))
        except Exception as e:
            print(e)
        finally:
            os.environ["PATH"] = str(orig_PATH)
            os.environ["GDAL_DATA"] = str(orig_GDAL_DATA)

    def create_mcu_model(self, dem):
        print('creating mcu model...')
        extents_gdf = gpd.read_file(self.dirs['dem_extents_path'])
        extents_gdf.crs = self.UTM
        clipped_vdatum_gdf = gpd.overlay(extents_gdf, self.vdatum_gdf, 
                                         how='intersection')
        region_geom = list(clipped_vdatum_gdf.geometry)
        region_mcu = clipped_vdatum_gdf.MCU_cm / 100
        mcu_model = features.rasterize(zip(region_geom, region_mcu), 
                                       out_shape=dem.shape, fill=np.nan, 
                                       transform=dem.transform)
        with rasterio.open(
                self.dirs['mcu_path'], 'w',
                driver='GTiff', 
                dtype=rasterio.float64, 
                count=1, 
                width=dem.width, 
                height=dem.height, 
                transform=dem.transform,
                nodata=np.nan) as dst:
            dst.write(mcu_model, indexes=1)

    def get_null_path(self):
        for i in range(self.num_chunks): 
            start_ind = i * self.chunk_size
            end_ind = start_ind + self.chunk_size
            xyz_chunk = self.xyz[start_ind:end_ind, :]
            null_path = self.dirs['null_ascii_path'].parent / self.dirs['null_ascii_path'].name.replace('.txt', f'{i}.txt')
            pd.DataFrame(xyz_chunk).to_csv(null_path, index=False)
            yield null_path

    def create_sep_model_multiprocess(self):
        p = pp.ProcessPool(4)
        for _ in tqdm(p.imap(self.create_sep_model, self.get_null_path()), 
                      total=self.num_chunks, ascii=True):
            pass
        p.close()
        p.join()


def get_dirs(dirs_path):
    with open(dirs_path) as json_file: 
        dirs = json.load(json_file) 

    paths = {}
    for k, v in dirs.items():
        paths[k] = Path(v)

    return paths


def set_env_vars(env_name):
    user_dir = os.path.expanduser('~')
    conda_dir = Path(user_dir).joinpath('AppData', 'Local', 
                                        'Continuum', 'anaconda3')
    env_dir = conda_dir / 'envs' / env_name
    share_dir = env_dir / 'Library' / 'share'
    script_path = conda_dir / 'Scripts'
    gdal_data_path = share_dir / 'gdal'
    proj_lib_path = share_dir

    if script_path.name not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep + str(script_path)
    os.environ["GDAL_DATA"] = str(gdal_data_path)
    os.environ["PROJ_LIB"] = str(proj_lib_path)


if __name__ == '__main__':
        set_env_vars('lidog_v2')
        py_dir = Path(os.path.dirname(os.path.realpath(__file__)))
        dirs = get_dirs(py_dir / 'dirs.json')

        jar_dir_parts = ('lib', 'java_home', 'openjdk-11.0.2', 'bin', 'java')
        dirs['java_dir'] = dirs['jar_dir'].joinpath(*jar_dir_parts)

        params = {
            'ihorz': 'NAD83_2011:utm:m:17',
            'ivert': 'NAD83_2011:m:height',
            'ohorz': 'NAD83_2011:utm:m:17',
            'overt': 'mllw:m:sounding:GEOID12B',
        }

        dem = rasterio.open(dirs['dem_path'])

        v_trans = VerticalTransform(dirs, params)
        v_trans.create_mcu_model(dem)
        v_trans.create_null_ascii(dem)
        v_trans.create_sep_model_multiprocess()
        v_trans.merge_sep_models(dem)
