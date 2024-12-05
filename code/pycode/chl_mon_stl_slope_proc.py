
from scipy.stats import linregress
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

ras_downsampled = ras[:, ::50, ::50]  # Select every 10th value for both latitude and longitude

# Initialize arrays to hold the slopes and p-values for the seascli component
seascli_slopes = np.zeros((ras_downsampled.shape[1], ras_downsampled.shape[2]))
seascli_p_values = np.ones((ras_downsampled.shape[1], ras_downsampled.shape[2]))

# Perform linear regression for the seascli component of each (lat, lon) cell
for i in range(ras.shape[1]):  # Loop over rows (latitude)
    for j in range(ras.shape[2]):  # Loop over columns (longitude)
        # Extract time series for the (i, j) cell
        ts = ras[:, i, j]
        
        # Skip cells with all NaN values
        if np.all(np.isnan(ts)):
            seascli_slopes[i, j] = np.nan
            seascli_p_values[i, j] = np.nan
            continue
        
        # Convert to pandas Series
        pix = pd.Series(ts, index=pd.date_range("1-1-1998", periods=len(ts), freq="ME"), name="chl")

        # Apply STL decomposition
        stl = STL(pix.ffill(), seasonal=13, robust=True)
        fit = stl.fit()
        
        # Compute the seascli component: trend + residual
        # trend = fit.trend #.dropna()  # Extract trend component
        resid = fit.resid #.dropna()  # Extract residual (anomaly) component
        
        # Not sure all of this is necissary since now using fffill
        # Ensure both components align in time
        # common_index = trend.index.intersection(residual.index)
        # seascli = trend.loc[common_index] + residual.loc[common_index]
        # seascli = trend + residual
        # Perform linear regression on the seascli component
        # x = (seascli.index - seascli.index[0]).days  # Time in days as the independent variable
        x = (resid.index - resid.index[0]).days  # Time in days as the independent variable
        # y = seascli.values
        y = resid.values
        
        # Compute slope and p-value using linregress
        slope, _, _, p_value, _ = linregress(x, y)
        seascli_slopes[i, j] = slope  # Store the slope
        seascli_p_values[i, j] = p_value  # Store the p-value

# The `seascli_slopes` and `seascli_p_values` arrays contain the regression results for the seascli component

# Create a new NetCDF file
ds = nc.Dataset('chla_stl_slope_20240902.nc', 'w', format='NETCDF4')

# Create dimensions
ds.createDimension('time', len(time))
ds.createDimension('latitude', len(lat))
ds.createDimension('longitude', len(lon))

# Create variables
times = ds.createVariable('time', 'f4', ('time',))
latitudes_var = ds.createVariable('latitude', 'f4', ('latitude',))
longitudes_var = ds.createVariable('longitude', 'f4', ('longitude',))
data_slope = ds.createVariable('slope', 'f4', ('latitude', 'longitude',), fill_value=np.nan)
data_pval = ds.createVariable('pval', 'f4', ('latitude', 'longitude',), fill_value=np.nan)

# Assign data to variables
times[:] = time
latitudes_var[:] = lat
longitudes_var[:] = lon
data_slope[:, :, :] = seascli_slopes
data_pval[:, :, :] = seascli_p_values

# Optionally, add attributes
times.units = 'days since 1900-01-01'
latitudes_var.units = 'degrees north'
longitudes_var.units = 'degrees east'
data_pval.units = 'fake'  
data_pval.long_name = 'significance of fit to resid plus trend.'
data_slope.units = 'mg m3 per yr'  
data_slope.long_name = 'slope of resid plus trend.'
# Close the file
ds.close()

print("NetCDF file created successfully.")