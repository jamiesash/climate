import numpy as np
from netCDF4 import Dataset

# -- Load data --
file_id = Dataset('/home/jamesash/koa_scratch/chla_day_med_month_1997_2025_20260328.nc')
ras  = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
lat  = file_id.variables["latitude"][:].copy()
lon  = file_id.variables["longitude"][:].copy()
time = file_id.variables["time"][:].copy()
time_units    = file_id.variables["time"].units
time_calendar = getattr(file_id.variables["time"], 'calendar', 'standard')
file_id.close()

ras[ras > 50] = np.nan

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# -- Keep only summer months: June (6) through October (10) --
months = date_vector.astype('datetime64[M]').astype(int) % 12 + 1
summer_mask = (months >= 6) & (months <= 10)

ras_summer  = ras[summer_mask, :, :]
time_summer = time[summer_mask]
date_summer = date_vector[summer_mask]

print(f'Total time steps: {len(date_vector)}')
print(f'Summer time steps: {len(date_summer)}')
print(f'First: {str(date_summer[0])[:10]}, Last: {str(date_summer[-1])[:10]}')

# -- Save summer-only data as NetCDF --
out_path = '/home/jamesash/koa_scratch/chla_day_med_month_summer_20260328.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

out.createDimension('time', ras_summer.shape[0])
out.createDimension('latitude', len(lat))
out.createDimension('longitude', len(lon))

time_var = out.createVariable('time', 'f8', ('time',))
time_var[:] = time_summer
time_var.units = time_units
time_var.calendar = time_calendar

lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

anom_var = out.createVariable('CHL_anom', 'f4', ('time', 'latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
anom_var[:] = ras_summer
anom_var.units = 'mg m^-3'
anom_var.long_name = 'Chlorophyll anomaly (daily median monthly removed), summer only (Jun-Oct)'

out.description = 'CHL anomaly (daily median monthly), summer months only (Jun-Oct)'
out.source_file = 'chla_day_med_month_1997_2025_20260328.nc'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')