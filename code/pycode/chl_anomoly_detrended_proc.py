from netCDF4 import Dataset
import numpy as np
from scipy.ndimage import uniform_filter1d

# -- Load data --
file_id = Dataset('../../data/chl/chl_1997_2025_day_l3_20260327.nc')
# ras  = file_id.variables["CHL"][:]
ras  = file_id.variables["CHL"][:].filled(np.nan).astype('float64')
lat  = file_id.variables["latitude"][:]
lon  = file_id.variables["longitude"][:]
time = file_id.variables["time"][:]
time_units = file_id.variables["time"].units
time_calendar = getattr(file_id.variables["time"], 'calendar', 'standard')
file_id.close()

# ras[ras > 100] = np.nan

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# ============================================================================
# Step 1: Remove seasonal cycle (monthly median climatology)
# ============================================================================
months = date_vector.astype('datetime64[M]').astype(int) % 12 + 1

monthly_clim = np.zeros((12, ras.shape[1], ras.shape[2]))
for m in range(1, 13):
    monthly_clim[m - 1, :, :] = np.nanmedian(ras[months == m, :, :], axis=0)

ras_anom = np.zeros_like(ras)
for i in range(ras.shape[0]):
    ras_anom[i, :, :] = ras[i, :, :] - monthly_clim[months[i] - 1, :, :]

# ============================================================================
# Step 2: Remove long-term trend (running mean)
# ============================================================================
# Window size in days — 1825 gives a 5-year running mean
# This smooths out interannual variability, leaving only the long-term trend
window = 1825

# Replace NaN with 0 for filtering, track valid counts
ras_filled = np.where(np.isnan(ras_anom), 0, ras_anom)
valid_count = np.where(np.isnan(ras_anom), 0, 1).astype('float64')

# Running mean along time axis (axis=0), handles edges with reflect
trend = uniform_filter1d(ras_filled, size=window, axis=0, mode='reflect')
count = uniform_filter1d(valid_count, size=window, axis=0, mode='reflect')

# Normalize by valid count to get true mean (correcting for NaN gaps)
count[count == 0] = np.nan
trend = trend / count

# Subtract the long-term trend
ras_detrended = ras_anom - trend

# Restore NaN where original was NaN
ras_detrended[np.isnan(ras_anom)] = np.nan

print(f'Anomaly range: {np.nanmin(ras_anom):.4f} to {np.nanmax(ras_anom):.4f}')
print(f'Trend range: {np.nanmin(trend):.4f} to {np.nanmax(trend):.4f}')
print(f'Detrended range: {np.nanmin(ras_detrended):.4f} to {np.nanmax(ras_detrended):.4f}')

# ============================================================================
# Save
# ============================================================================
out_path = '/home/jamesash/koa_scratch/chla_day_deseason_detrend_20260422.nc'
out = Dataset(out_path, 'w', format='NETCDF4')

out.createDimension('time', ras_detrended.shape[0])
out.createDimension('latitude', len(lat))
out.createDimension('longitude', len(lon))
out.createDimension('month', 12)

time_var = out.createVariable('time', 'f8', ('time',))
time_var[:] = time[:]
time_var.units = time_units
time_var.calendar = time_calendar

lat_var = out.createVariable('latitude', 'f4', ('latitude',))
lat_var[:] = lat
lat_var.units = 'degrees_north'

lon_var = out.createVariable('longitude', 'f4', ('longitude',))
lon_var[:] = lon
lon_var.units = 'degrees_east'

# Detrended anomaly
detrend_var = out.createVariable('CHL_anom', 'f4', ('time', 'latitude', 'longitude'),
                                 fill_value=np.nan, zlib=True)
detrend_var[:] = ras_detrended
detrend_var.units = 'mg m^-3'
detrend_var.long_name = 'Chlorophyll anomaly (seasonal and long-term trend removed)'

# Long-term trend (for inspection)
trend_var = out.createVariable('CHL_trend', 'f4', ('time', 'latitude', 'longitude'),
                               fill_value=np.nan, zlib=True)
trend_var[:] = trend
trend_var.units = 'mg m^-3'
trend_var.long_name = f'Long-term trend ({window}-day running mean of seasonal anomaly)'

# Monthly climatology
clim_var = out.createVariable('CHL_clim', 'f4', ('month', 'latitude', 'longitude'),
                              fill_value=np.nan, zlib=True)
clim_var[:] = monthly_clim
clim_var.units = 'mg m^-3'
clim_var.long_name = 'Monthly median climatology (12 months)'

out.description = 'CHL with seasonal cycle and long-term trend removed'
out.source_file = 'chl_1997_2025_day_l3_20260327.nc'
out.deseason_method = 'Monthly median climatology subtracted'
out.detrend_method = f'{window}-day running mean subtracted from seasonal anomaly'
out.history = f'Created {np.datetime64("today")}'

out.close()
print(f'Saved: {out_path}')
