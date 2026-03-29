import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import geopandas as gpd

# -- Load data --
file_id = Dataset('../../data/chl/chl_1997_2025_day_l3_20260327.nc')
ras  = file_id.variables["CHL"][:]
lat  = file_id.variables["latitude"][:]
lon  = file_id.variables["longitude"][:]
time = file_id.variables["time"][:]
time_units = file_id.variables["time"].units       # preserve original units
time_calendar = getattr(file_id.variables["time"], 'calendar', 'standard')
file_id.close()

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# -- Extract month for each time step (works for daily or monthly) --
months = date_vector.astype('datetime64[M]').astype(int) % 12 + 1

# -- Compute median monthly climatology (12 months, each pixel) --
monthly_clim = np.zeros((12, ras.shape[1], ras.shape[2]))
for m in range(1, 13):
    monthly_clim[m - 1, :, :] = np.nanmedian(ras[months == m, :, :], axis=0)

# -- Subtract the climatology from each time step --
ras_anom = np.zeros_like(ras)
for i in range(ras.shape[0]):
    ras_anom[i, :, :] = ras[i, :, :] - monthly_clim[months[i] - 1, :, :]

# -- Save anomaly raster as NetCDF --
out_path = '/home/jamesash/koa_scratch/chla_day_med_month_1997_2025_20260328.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

# Dimensions
out.createDimension('time', ras_anom.shape[0])
out.createDimension('latitude', len(lat))
out.createDimension('longitude', len(lon))
out.createDimension('month', 12)

# Time variable — preserve original units and values
time_var = out.createVariable('time', 'f8', ('time',))
time_var[:] = time[:]
time_var.units = time_units
time_var.calendar = time_calendar

# Coordinate variables
lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

# Anomaly data
anom_var = out.createVariable('CHL_anom', 'f4', ('time', 'latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
anom_var[:] = ras_anom
anom_var.units = 'mg m^-3'
anom_var.long_name = 'Chlorophyll anomaly (observed minus monthly median climatology)'

# Monthly climatology
clim_var = out.createVariable('CHL_clim', 'f4', ('month', 'latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
clim_var[:] = monthly_clim
clim_var.units = 'mg m^-3'
clim_var.long_name = 'Monthly median climatology (12 months)'

# Global attributes
out.description = 'CHL monthly median anomaly and climatology'
out.source_file = 'chl_1998_2025_l4_month_multi_4k.nc'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')
