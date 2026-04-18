import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.ndimage import center_of_mass

# -- Load cropped and masked anomaly data --
file_id = Dataset('/home/jamesash/koa_scratch/chl_anomaly_cropped_masked_20260328.nc')
# ras  = file_id.variables["CHL_anom"][:].copy()
ras  = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
lat  = file_id.variables["latitude"][:]
lon  = file_id.variables["longitude"][:]
time = file_id.variables["time"][:]
file_id.close()

# fix the 999 fill value
ras[ras > 50] = np.nan

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# -- Compute per-pixel median and MAD, then mask extremes --
pixel_median = np.nanmedian(ras, axis=0)
pixel_mad    = np.nanmedian(np.abs(ras - pixel_median), axis=0)
# Scaled MAD (consistent estimator of std for normal data)
pixel_mad_scaled = pixel_mad * 1.4826
ras_mask = np.zeros_like(ras)
ras_mask[ras > pixel_median + pixel_mad_scaled] = 1

ras_extreme = np.where(ras_mask == 1, ras, np.nan)

# -- Define bloom date windows --
start_dates = [
    '1998-07-08', '1999-08-01', '2000-08-15', '2001-07-18',
    '2002-06-15', '2003-07-26', '2004-06-03', '2005-07-02', '2006-07-22',
    '2007-08-15', '2008-08-10', '2009-06-13', '2010-07-01', '2011-09-07',
    '2012-09-01', '2013-07-10', '2014-06-14', '2015-08-01', '2016-08-02',
    '2017-07-15', '2018-07-16', '2019-06-11', '2020-08-10', '2021-08-01',
    '2022-07-01', '2023-09-21', '2024-06-27', '2025-07-22'
]
end_dates = [
    '1998-10-06', '1999-10-12', '2000-11-15', '2001-08-04',
    '2002-09-10', '2003-10-08', '2004-08-11', '2005-10-02', '2006-11-11',
    '2007-11-07', '2008-10-22', '2009-10-10', '2010-09-28', '2011-11-12',
    '2012-11-06', '2013-10-02', '2014-09-01', '2015-10-01', '2016-09-25',
    '2017-09-05', '2018-10-27', '2019-09-08', '2020-09-08', '2021-10-01',
    '2022-10-10', '2023-10-28', '2024-09-19', '2025-10-29'
]

start_dates = np.array(start_dates, dtype='datetime64[D]')
end_dates   = np.array(end_dates, dtype='datetime64[D]')

# -- Compute magnitude and center of mass for each bloom --
results = []
for s, e in zip(start_dates, end_dates):
    tmask = (date_vector >= s) & (date_vector <= e)
    subset_dates = date_vector[tmask]
    subset_extreme = ras_extreme[tmask, :, :]

    # Magnitude
    # u = np.nanmean(subset_extreme)
    # o = np.nanstd(subset_extreme)
    # mag = u + 2 * o
    mag = np.nanpercentile(subset_extreme, 95)

    # Center of mass
    subset_com = np.where(np.isnan(subset_extreme), 0, subset_extreme)

    # power amplify to give higher values more wieght. 
    power = 2
    subset_weighted = subset_com ** power

    if np.sum(subset_com) == 0:
        results.append({
            'start': str(s),
            'end': str(e),
            'center_date': 'N/A',
            'center_lat': np.nan,
            'center_lon': np.nan,
            'magnitude': np.nan
        })
        continue

    t_idx, lat_idx, lon_idx = center_of_mass(subset_weighted)
    t_round = min(int(round(t_idx)), len(subset_dates) - 1)
    center_date = str(subset_dates[t_round])[:10]

    results.append({
        'start': str(s),
        'end': str(e),
        'center_date': center_date,
        'center_lat': float(lat[int(round(lat_idx))]),
        'center_lon': float(lon[int(round(lon_idx))]),
        'magnitude': mag
    })

# -- Save to CSV --
df = pd.DataFrame(results)
out_csv = '/home/jamesash/koa_scratch/bloom_summary_20260328.csv'
df.to_csv(out_csv, index=False)
