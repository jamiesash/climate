import numpy as np
from netCDF4 import Dataset

# -- Load data --
file_id = Dataset('/home/jamesash/koa_scratch/chla_day_med_month_1997_2025_20260328.nc')
ras  = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
lat  = file_id.variables["latitude"][:].copy()
lon  = file_id.variables["longitude"][:].copy()
time = file_id.variables["time"][:].copy()
file_id.close()

ras[ras > 50] = np.nan

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# -- Keep only summer months: June (6) through October (10) --
months = date_vector.astype('datetime64[M]').astype(int) % 12 + 1
summer_mask = (months >= 6) & (months <= 10)

ras_summer = ras[summer_mask, :, :]
date_summer = date_vector[summer_mask]

print(f'Total time steps: {len(date_vector)}')
print(f'Summer time steps: {len(date_summer)}')
print(f'First: {str(date_summer[0])[:10]}, Last: {str(date_summer[-1])[:10]}')

# -- Compute mean per grid cell across all summer time steps --
ras_mean = np.nanmean(ras_summer, axis=0)

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
mean_var.long_name = 'Mean summer (Jun-Oct) chlorophyll anomaly across all years'

out.description = 'Temporal mean of summer CHL anomaly (Jun-Oct) per grid cell'
out.source_file = 'chla_day_med_month_1997_2025_20260328.nc'
out.summer_months = 'June through October (6-10)'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')