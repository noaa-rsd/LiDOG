set /P env_dir="Enter conda environments directory: "
call %UserProfile%\AppData\Local\Continuum\anaconda3\condabin\conda.bat create --clone "C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3" --prefix %env_dir%\lidog
call %UserProfile%\AppData\Local\Continuum\anaconda3\condabin\conda.bat activate lidog
call %UserProfile%\AppData\Local\Continuum\anaconda3\condabin\conda.bat install -c conda-forge rasterio shapely pyproj