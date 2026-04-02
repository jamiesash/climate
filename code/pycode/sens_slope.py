import numpy as np
from netCDF4 import Dataset
from scipy.stats import mstats

# -- Load summer anomaly data --
file_id = Dataset('/home/jamesash/koa_scratch/chl_anomaly_summer_20260328.nc')
ras  = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
lat  = file_id.variables["latitude"][:].copy()
lon  = file_id.variables["longitude"][:].copy()
time = file_id.variables["time"][:].copy()
time_units    = file_id.variables["time"].units
time_calendar = getattr(file_id.variables["time"], 'calendar', 'standard')
file_id.close()

ras[ras > 100] = np.nan

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# Time as numeric values (days since first observation)
time_numeric = (date_vector - date_vector[0]).astype('timedelta64[D]').astype(float)

# -- Sen's slope at each pixel --
nlat = len(lat)
nlon = len(lon)
slope_map     = np.full((nlat, nlon), np.nan)
intercept_map = np.full((nlat, nlon), np.nan)
lo_slope_map  = np.full((nlat, nlon), np.nan)
hi_slope_map  = np.full((nlat, nlon), np.nan)

total_pixels = nlat * nlon
print(f'Running Sen slope on {total_pixels} pixels...')

for i in range(nlat):
    if i % 50 == 0:
        print(f'  Row {i}/{nlat}')
    for j in range(nlon):
        ts = ras[:, i, j]
        valid = ~np.isnan(ts)

        if np.sum(valid) < 10:
            continue

        slope, intercept, lo, hi = mstats.theilslopes(ts[valid], time_numeric[valid])
        slope_map[i, j]     = slope
        intercept_map[i, j] = intercept
        lo_slope_map[i, j]  = lo
        hi_slope_map[i, j]  = hi

print('Done.')

# -- Save results as NetCDF --
out_path = '/home/jamesash/koa_scratch/chl_senslope_summer_20260328.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

out.createDimension('latitude', nlat)
out.createDimension('longitude', nlon)

lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

slope_var = out.createVariable('slope', 'f4', ('latitude', 'longitude'),
                               fill_value=np.nan, zlib=True)
slope_var[:] = slope_map
slope_var.units = 'mg m^-3 per day'
slope_var.long_name = 'Sen slope (CHL anomaly vs time)'

intercept_var = out.createVariable('intercept', 'f4', ('latitude', 'longitude'),
                                   fill_value=np.nan, zlib=True)
intercept_var[:] = intercept_map
intercept_var.units = 'mg m^-3'
intercept_var.long_name = 'Sen slope intercept'

lo_var = out.createVariable('slope_lo', 'f4', ('latitude', 'longitude'),
                            fill_value=np.nan, zlib=True)
lo_var[:] = lo_slope_map
lo_var.units = 'mg m^-3 per day'
lo_var.long_name = 'Sen slope lower 95% confidence bound'

hi_var = out.createVariable('slope_hi', 'f4', ('latitude', 'longitude'),
                            fill_value=np.nan, zlib=True)
hi_var[:] = hi_slope_map
hi_var.units = 'mg m^-3 per day'
hi_var.long_name = 'Sen slope upper 95% confidence bound'

out.description = 'Sen slope trend on summer CHL anomaly (Jun-Oct)'
out.source_file = 'chl_anomaly_summer_20260328.nc'
out.time_units = 'slope is per day; multiply by 365.25 for per year'
out.min_valid_obs = '10'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')