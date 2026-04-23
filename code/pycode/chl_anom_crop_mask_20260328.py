import numpy as np
import geopandas as gpd
from netCDF4 import Dataset
from shapely.ops import polygonize, unary_union
from shapely.geometry import Polygon, MultiPolygon
from matplotlib.path import Path

# -- Load anomaly data --
file_id = Dataset('/home/jamesash/koa_scratch/chla_day_deseason_detrend_20260423.nc')
ras      = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
clim     = file_id.variables["CHL_clim"][:]
lat      = file_id.variables["latitude"][:]
lon      = file_id.variables["longitude"][:]
time     = file_id.variables["time"][:]
time_units    = file_id.variables["time"].units
time_calendar = getattr(file_id.variables["time"], 'calendar', 'standard')
file_id.close()

# -- Crop to region of interest: 18N-35N, 170W-130W --
lat_mask = (lat >= 18) & (lat <= 35)
lon_mask = (lon >= -170) & (lon <= -130)
lat  = lat[lat_mask]
lon  = lon[lon_mask]
ras  = ras[:, lat_mask, :][:, :, lon_mask]
clim = clim[:, lat_mask, :][:, :, lon_mask]

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# -- Load and build the Hawaii CZ polygons --
gdf = gpd.read_file('../../data/eez/USMaritimeLimitsNBoundaries.shp')
hawaii = gdf[gdf['REGION'] == 'Hawaiian Islands']
cz_lines = hawaii[hawaii['CZ'] == 1.0]
merged_lines = unary_union(cz_lines.geometry)
cz_polys = list(polygonize(merged_lines))
cz_combined = unary_union(cz_polys)

# -- Create a 2D mask on the lat/lon grid --
lon_grid, lat_grid = np.meshgrid(lon, lat)
points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))

if isinstance(cz_combined, Polygon):
    cz_combined = MultiPolygon([cz_combined])

inside = np.zeros(len(points), dtype=bool)
for poly in cz_combined.geoms:
    p = Path(np.array(poly.exterior.coords))
    inside |= p.contains_points(points)
inside = inside.reshape(lon_grid.shape)

# -- Apply mask: set pixels inside the CZ to NaN --
ras[:, inside] = np.nan

# -- Save cropped and masked anomaly as NetCDF --
out_path = '/home/jamesash/koa_scratch/chl_anomaly_cropped_masked_20260423.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

# Dimensions
out.createDimension('time', ras.shape[0])
out.createDimension('latitude', len(lat))
out.createDimension('longitude', len(lon))
out.createDimension('month', 12)

# Time
time_var = out.createVariable('time', 'f8', ('time',))
time_var[:] = time[:]
time_var.units = time_units
time_var.calendar = time_calendar

# Coordinates
lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

# Anomaly
anom_var = out.createVariable('CHL_anom', 'f4', ('time', 'latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
anom_var[:] = ras
anom_var.units = 'mg m^-3'
anom_var.long_name = 'Chlorophyll anomaly, cropped and CZ masked'

# Climatology (cropped)
clim_var = out.createVariable('CHL_clim', 'f4', ('month', 'latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
clim_var[:] = clim
clim_var.units = 'mg m^-3'
clim_var.long_name = 'Monthly median climatology (12 months), cropped'

# CZ mask
mask_var = out.createVariable('CZ_mask', 'i1', ('latitude', 'longitude'), zlib=True)
mask_var[:] = inside.astype(np.int8)
mask_var.long_name = 'Contiguous zone mask (1 = inside CZ, 0 = open ocean)'

# Global attributes
out.description = 'CHL anomaly cropped to 18-35N, 170-130W with Hawaii CZ masked'
out.source_file = 'chl_anomaly_20260328.nc'
out.spatial_bounds = '18N-35N, 170W-130W'
out.mask_applied = 'Hawaii Contiguous Zone (CZ) set to NaN'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')
