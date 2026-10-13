"""
Download and quality-control Louisiana precipitation records from NOAA GHCN-Daily.

Reads the fixed-width .dly files directly, keeps PRCP, screens on NOAA's own
quality flags, aggregates to monthly totals, and refuses any month that is not
sufficiently complete. Writes data/monthly_precip.csv.

    python fetch_ghcn.py
"""
import io, os, urllib.request
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
RAW = os.path.join(DATA, 'raw')
os.makedirs(RAW, exist_ok=True)

BASE = 'https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/%s.dly'

# Eight long-record stations spread across Louisiana's climate regions.
STATIONS = {
    'USC00167344': 'Plain Dealing (NW)',
    'USC00164355': 'Homer (N)',
    'USC00169806': 'Winnsboro (NE)',
    'USC00160098': 'Alexandria (Central)',
    'USW00013976': 'Lafayette Regional AP (SW)',
    'USC00166657': 'New Iberia (S Central)',
    'USC00162151': 'Covington (SE)',
    'USW00012930': 'New Orleans Audubon (SE Coastal)',
}

START_YEAR = 1950          # common period all eight stations cover
MIN_DAYS_FRAC = 0.90       # a month needs this share of valid days to be used

# GHCN-Daily quality flags. Any non-blank value here means the observation
# failed one of NOAA's checks, so it is discarded rather than treated as zero.
BAD_QFLAGS = set('DGIKLMNORSTWXZ')


def station_meta():
    """Latitude, longitude and elevation from ghcnd-stations.txt."""
    meta = {}
    path = os.path.join(DATA, 'ghcnd-stations.txt')
    for line in io.open(path, encoding='utf-8'):
        sid = line[0:11]
        if sid in STATIONS:
            meta[sid] = {
                'lat': float(line[12:20]), 'lon': float(line[21:30]),
                'elev_m': float(line[31:37]), 'ghcn_name': line[41:71].strip(),
            }
    return meta


def parse_dly(path, sid):
    """Yield (year, month, day, prcp_mm) for valid PRCP observations."""
    for line in io.open(path, encoding='utf-8'):
        if line[17:21] != 'PRCP':
            continue
        year, month = int(line[11:15]), int(line[15:17])
        if year < START_YEAR:
            continue
        for d in range(31):
            off = 21 + d * 8
            raw = line[off:off + 5].strip()
            qflag = line[off + 6:off + 7]
            if raw == '' or raw == '-9999':
                continue
            if qflag.strip() and qflag in BAD_QFLAGS:
                continue                       # failed a NOAA quality check
            yield year, month, d + 1, int(raw) / 10.0   # tenths of mm -> mm


rows, qc = [], []
meta = station_meta()

for sid, label in STATIONS.items():
    local = os.path.join(RAW, sid + '.dly')
    if not os.path.exists(local):
        print('downloading %s  %s' % (sid, label))
        urllib.request.urlretrieve(BASE % sid, local)
    obs = pd.DataFrame(parse_dly(local, sid), columns=['year', 'month', 'day', 'prcp_mm'])
    if obs.empty:
        print('  !! no usable PRCP for %s' % sid)
        continue

    g = obs.groupby(['year', 'month'])['prcp_mm']
    monthly = g.agg(prcp_mm='sum', n_days='size').reset_index()
    days_in_month = pd.to_datetime(
        dict(year=monthly.year, month=monthly.month, day=1)).dt.days_in_month
    monthly['completeness'] = monthly['n_days'] / days_in_month

    kept = monthly[monthly['completeness'] >= MIN_DAYS_FRAC].copy()
    dropped = len(monthly) - len(kept)

    kept['station_id'] = sid
    kept['station'] = label
    for k, v in meta.get(sid, {}).items():
        kept[k] = v
    rows.append(kept)

    qc.append({
        'station_id': sid, 'station': label,
        'first_year': int(monthly.year.min()), 'last_year': int(monthly.year.max()),
        'months_total': len(monthly), 'months_kept': len(kept),
        'months_dropped_incomplete': dropped,
        'pct_kept': round(100 * len(kept) / len(monthly), 1),
        'daily_obs_used': int(obs.shape[0]),
    })
    print('  %s: %d of %d months pass completeness (%.1f%%)'
          % (label, len(kept), len(monthly), 100 * len(kept) / len(monthly)))

out = pd.concat(rows, ignore_index=True)
out = out[['station_id', 'station', 'lat', 'lon', 'elev_m',
           'year', 'month', 'prcp_mm', 'n_days', 'completeness']]
out = out.sort_values(['station_id', 'year', 'month']).reset_index(drop=True)
out.to_csv(os.path.join(DATA, 'monthly_precip.csv'), index=False)

qcdf = pd.DataFrame(qc)
qcdf.to_csv(os.path.join(DATA, 'station_qc_summary.csv'), index=False)

print()
print(qcdf.to_string(index=False))
print()
print('monthly records written : %d' % len(out))
print('period                  : %d-%d' % (out.year.min(), out.year.max()))
print('stations                : %d' % out.station_id.nunique())
