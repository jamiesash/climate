import numpy as np
from netCDF4 import Dataset
from scipy.stats import kendalltau

# -- Load summer anomaly data --
file_id = Dataset('/home/jamesash/koa_scratch/chla_day_summer_20260403.nc')
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

# Time as numeric values (days since base date) for correlation
time_numeric = (date_vector - date_vector[0]).astype('timedelta64[D]').astype(float)

# -- Kendall's Tau at each pixel --
nlat = len(lat)
nlon = len(lon)
tau_map = np.full((nlat, nlon), np.nan)
pval_map = np.full((nlat, nlon), np.nan)

total_pixels = nlat * nlon
print(f'Running Kendall Tau on {total_pixels} pixels...')

for i in range(nlat):
    if i % 50 == 0:
        print(f'  Row {i}/{nlat}')
    for j in range(nlon):
        ts = ras[:, i, j]
        valid = ~np.isnan(ts)

        # Need at least 10 valid observations
        if np.sum(valid) < 10:
            continue

        tau, pval = kendalltau(time_numeric[valid], ts[valid])
        tau_map[i, j] = tau
        pval_map[i, j] = pval

print('Done.')

# -- Save results as NetCDF --
out_path = '/home/jamesash/koa_scratch/chl_kendall_20260405.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

out.createDimension('latitude', nlat)
out.createDimension('longitude', nlon)

lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

tau_var = out.createVariable('tau', 'f4', ('latitude', 'longitude'),
                             fill_value=np.nan, zlib=True)
tau_var[:] = tau_map
tau_var.long_name = 'Kendall Tau correlation coefficient (CHL anomaly vs time)'

pval_var = out.createVariable('pvalue', 'f4', ('latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
pval_var[:] = pval_map
pval_var.long_name = 'Kendall Tau p-value (two-sided)'

out.description = 'Kendall Tau trend test on summer CHL anomaly (Jun-Oct)'
out.source_file = 'chl_anomaly_summer_20260328.nc'
out.min_valid_obs = '10'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')
