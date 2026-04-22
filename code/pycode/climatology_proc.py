import numpy as np
from netCDF4 import Dataset

# -- Load summer anomaly data --
file_id = Dataset('/home/jamesash/koa_scratch/chla_day_summer_20260403.nc')
ras  = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
lat  = file_id.variables["latitude"][:].copy()
lon  = file_id.variables["longitude"][:].copy()
file_id.close()

ras[ras > 50] = np.nan

# -- Compute mean per grid cell across all time steps --
ras_mean = np.nanmean(ras, axis=0)  # shape: (lat, lon)

print(f'Input shape: {ras.shape}')
print(f'Output shape: {ras_mean.shape}')
print(f'Value range: {np.nanmin(ras_mean):.4f} to {np.nanmax(ras_mean):.4f}')

# -- Save as NetCDF --
out_path = '/home/jamesash/koa_scratch/chla_day_summer_mean_20260403.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

out.createDimension('latitude', len(lat))
out.createDimension('longitude', len(lon))

lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

mean_var = out.createVariable('CHL_anom_mean', 'f4', ('latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
mean_var[:] = ras_mean
mean_var.units = 'mg m^-3'
mean_var.long_name = 'Mean summer chlorophyll anomaly across all years'

out.description = 'Temporal mean of summer CHL anomaly per grid cell'
out.source_file = 'chla_day_summer_20260403.nc'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')