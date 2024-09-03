import numpy as np
import netCDF4 as nc
import pandas as pd
from statsmodels.tsa.seasonal import STL

# Global Ocean Colour (Copernicus-GlobColour), Bio-Geo-Chemical, L4 (monthly and interpolated) from Satellite Observations (Near Real Time)
file_id = nc.Dataset('/home/jamie/projects/climate/data/chl/chl_1998_2023_l4_month_multi_4k.nc')
ras = file_id.variables["CHL"][:]
time = file_id.variables["time"][:]
lat = file_id.variables["latitude"][:]
lon = file_id.variables["longitude"][:]
file_id.close()

# Initialize an array to hold the slopes
tmp = np.zeros((ras.shape[0], ras.shape[1], ras.shape[2]))  # shape (m, n)

# Perform linear regression for each (m, n) cell
for i in range(ras.shape[1]):  # Loop over rows
    for j in range(ras.shape[2]):  # Loop over columns
        # Perform linear regression for the (i, j) cell over time
        ts = ras[:, i, j]
        pix = pd.Series(ts, index=pd.date_range("1-1-1998", periods=len(ts), freq="M"), name="chl")
        stl = STL(pix, seasonal = 13, robust = True)
        fit = stl.fit()
        tmp[:, i, j] = fit.resid  # Store the slope in the 2D array

# Create a new NetCDF file
ds = nc.Dataset('chla_stl_mon_20240902.nc', 'w', format='NETCDF4')

# Create dimensions
ds.createDimension('time', len(time))
ds.createDimension('latitude', len(lat))
ds.createDimension('longitude', len(lon))

# Create variables
times = ds.createVariable('time', 'f4', ('time',))
latitudes_var = ds.createVariable('latitude', 'f4', ('latitude',))
longitudes_var = ds.createVariable('longitude', 'f4', ('longitude',))
data_var = ds.createVariable('CHL', 'f4', ('time', 'latitude', 'longitude',), fill_value=np.nan)

# Assign data to variables
times[:] = time
latitudes_var[:] = lat
longitudes_var[:] = lon
data_var[:, :, :] = tmp

# Optionally, add attributes
times.units = 'days since 1900-01-01'
latitudes_var.units = 'degrees north'
longitudes_var.units = 'degrees east'
data_var.units = 'mg m3'  
data_var.long_name = 'seasonal signal from stl on monthly chl.'

# Close the file
ds.close()

print("NetCDF file created successfully.")