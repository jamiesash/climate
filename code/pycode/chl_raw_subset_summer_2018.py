import numpy as np
from netCDF4 import Dataset

# -- Load cropped and masked raw CHL data --
file_id = Dataset('/home/jamesash/koa_scratch/chl_raw_cropped_masked_20260328.nc')
ras  = file_id.variables["CHL"][:].filled(np.nan).astype('float64')
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

# -- Subset to 2018 bloom period --
bloom_mask = (date_vector >= np.datetime64('2018-06-19')) & \
             (date_vector <= np.datetime64('2018-10-30'))

ras_2018  = ras[bloom_mask, :, :]
time_2018 = time[bloom_mask]
date_2018 = date_vector[bloom_mask]

print(f'2018 bloom raw CHL time steps: {len(date_2018)}')
print(f'First: {str(date_2018[0])[:10]}, Last: {str(date_2018[-1])[:10]}')

# -- Save --
out_path = '/home/jamesash/koa_scratch/chl_raw_2018_bloom_20260328.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

out.createDimension('time', ras_2018.shape[0])
out.createDimension('latitude', len(lat))
out.createDimension('longitude', len(lon))

time_var = out.createVariable('time', 'f8', ('time',))
time_var[:] = time_2018
time_var.units = time_units
time_var.calendar = time_calendar

lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

chl_var = out.createVariable('CHL', 'f4', ('time', 'latitude', 'longitude'),
                             fill_value=np.nan, zlib=True)
chl_var[:] = ras_2018
chl_var.units = 'mg m^-3'
chl_var.long_name = 'Raw chlorophyll, 2018 bloom period'

out.description = '2018 bloom raw CHL subset (2018-06-19 to 2018-10-30)'
out.source_file = 'chl_raw_cropped_masked_20260328.nc'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')