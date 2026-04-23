import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.ndimage import center_of_mass

# -- Load cropped and masked anomaly data --
file_id = Dataset('/home/jamesash/koa_scratch/chl_anomaly_cropped_masked_20260423.nc')
ras  = file_id.variables["CHL_anom"][:].filled(np.nan).astype('float64')
cz_mask = file_id.variables["CZ_mask"][:].astype(bool)  # True = inside CZ
lat  = file_id.variables["latitude"][:].copy()
lon  = file_id.variables["longitude"][:].copy()
time = file_id.variables["time"][:].copy()
file_id.close()

# -- Load cropped and masked raw CHL data --
file_id = Dataset('/home/jamesash/koa_scratch/chl_raw_cropped_masked_20260328.nc')
ras_raw = file_id.variables["CHL"][:].filled(np.nan).astype('float64')
file_id.close()

# Replace fill values
ras[ras > 80] = np.nan
ras_raw[ras_raw > 80] = np.nan

# -- Build date vector --
timedelta_vector = (time * np.timedelta64(1, 'D')).astype('timedelta64[ns]')
base_date   = np.datetime64('1900-01-01')
date_vector = base_date + timedelta_vector

# -- Compute per-pixel mean and std on anomaly, then create mask --
# pixel_mean = np.nanmean(ras, axis=0)
# pixel_std  = np.nanstd(ras, axis=0)
# ras_mask = np.zeros_like(ras)
# ras_mask[ras > pixel_mean + pixel_std] = 1

# -- Compute per-pixel median and MAD, then mask extremes --
pixel_median = np.nanmedian(ras, axis=0)
pixel_mad    = np.nanmedian(np.abs(ras - pixel_median), axis=0)
# Scaled MAD (consistent estimator of std for normal data)
pixel_mad_scaled = pixel_mad * 1.4826
ras_mask = np.zeros_like(ras)
ras_mask[ras > pixel_median + pixel_mad_scaled] = 1

# Apply the anomaly-derived mask to BOTH datasets
ras_extreme_anom = np.where(ras_mask == 1, ras, np.nan)
ras_extreme_raw  = np.where(ras_mask == 1, ras_raw, np.nan)

# -- Precompute pixel areas (km²) accounting for latitude --
dlat = np.abs(np.median(np.diff(lat)))
dlon = np.abs(np.median(np.diff(lon)))
km_per_deg_lat = 111.32
lat_rad = np.deg2rad(lat)
pixel_area_km2 = (dlat * km_per_deg_lat) * (dlon * km_per_deg_lat * np.cos(lat_rad))

# -- Precompute ocean pixel count (excluding CZ) for cloud cover --
ocean_pixels = np.sum(~cz_mask)  # total non-CZ pixels

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

# -- Compute magnitude, center of mass, max area, and cloud cover --
results = []
for s, e in zip(start_dates, end_dates):
    tmask = (date_vector >= s) & (date_vector <= e)
    subset_dates = date_vector[tmask]
    subset_extreme_anom = ras_extreme_anom[tmask, :, :]
    subset_extreme_raw  = ras_extreme_raw[tmask, :, :]
    subset_mask = ras_mask[tmask, :, :]
    subset_raw  = ras_raw[tmask, :, :]

    # Magnitude from RAW CHL at extreme pixels
    mag = np.nanpercentile(subset_extreme_raw, 95)

    # Center of mass from ANOMALY extreme values
    subset_com = np.where(np.isnan(subset_extreme_anom), 0, subset_extreme_anom)

    # the uncoimmented one is not necissarily what is in the results. 
    # subset_weighted = subset_com**2
    subset_weighted = (np.exp(subset_com) -1)**2

    if np.sum(subset_com) == 0:
        results.append({
            'start': str(s),
            'end': str(e),
            'center_date': 'N/A',
            'center_lat': np.nan,
            'center_lon': np.nan,
            'magnitude': np.nan,
            'max_area_km2': np.nan,
            'max_area_date': 'N/A',
            'cloud_pct': np.nan
        })
        continue

    t_idx, lat_idx, lon_idx = center_of_mass(subset_weighted)
    t_round = min(int(round(t_idx)), len(subset_dates) - 1)
    center_date = str(subset_dates[t_round])[:10]

    # Area for each day in the window
    daily_areas = np.array([
        np.sum(subset_mask[t, :, :] * pixel_area_km2[:, np.newaxis])
        for t in range(subset_mask.shape[0])
    ])

    # Day with the largest bloom area
    max_t = np.argmax(daily_areas)
    max_area = daily_areas[max_t]
    max_area_date = str(subset_dates[max_t])[:10]

    # Cloud cover on the max area day (NaN in ocean pixels = cloud)
    day_slice = subset_raw[max_t, :, :]
    nan_in_ocean = np.isnan(day_slice) & (~cz_mask)
    cloud_pct = (np.sum(nan_in_ocean) / ocean_pixels) * 100.0

    results.append({
        'start': str(s),
        'end': str(e),
        'center_date': center_date,
        'center_lat': float(lat[int(round(lat_idx))]),
        'center_lon': float(lon[int(round(lon_idx))]),
        'magnitude': mag,
        'max_area_km2': float(max_area),
        'max_area_date': max_area_date,
        'cloud_pct': round(float(cloud_pct), 2)
    })

# -- Save to CSV --
df = pd.DataFrame(results)
out_csv = '/home/jamesash/koa_scratch/bloom_summary_exp_2p_1mad_20260423.csv'
df.to_csv(out_csv, index=False)
print(f'Saved: {out_csv}')
print(df.to_string())
