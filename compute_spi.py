"""
Standardised Precipitation Index and drought characteristics for the Louisiana
station network.

Method follows McKee et al. (1993). For each station, timescale and calendar
month, a two-parameter gamma distribution is fitted to the accumulated
precipitation, zero-precipitation months are handled with a mixed distribution,
and the cumulative probability is mapped to the standard normal. SPI is
therefore in units of standard deviations: -1 is the 16th percentile of the
local record for that month, -2 roughly the 2nd.

Writes  data/spi.csv, data/drought_events.csv,
        data/summary_by_station.csv, data/summary_by_season.csv

    python compute_spi.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

SCALES = [3, 6, 12]          # months of accumulation
DROUGHT_THRESHOLD = -1.0     # SPI at or below this counts as drought
SEASONS = {12: 'Winter (DJF)', 1: 'Winter (DJF)', 2: 'Winter (DJF)',
           3: 'Spring (MAM)', 4: 'Spring (MAM)', 5: 'Spring (MAM)',
           6: 'Summer (JJA)', 7: 'Summer (JJA)', 8: 'Summer (JJA)',
           9: 'Autumn (SON)', 10: 'Autumn (SON)', 11: 'Autumn (SON)'}


def classify(spi):
    if pd.isna(spi):
        return None
    if spi <= -2.0:
        return 'Extreme'
    if spi <= -1.5:
        return 'Severe'
    if spi <= -1.0:
        return 'Moderate'
    return 'None'


def spi_series(acc, months):
    """SPI for one station and timescale, fitted per calendar month."""
    out = pd.Series(np.nan, index=acc.index)
    for m in range(1, 13):
        sel = (months == m) & acc.notna()
        x = acc[sel]
        if len(x) < 30:                       # too short to fit a distribution
            continue
        zero = (x <= 0).sum()
        q = zero / len(x)                     # probability of a dry accumulation
        pos = x[x > 0]
        if len(pos) < 20:
            continue
        shape, loc, scale = stats.gamma.fit(pos, floc=0)
        cdf = q + (1 - q) * stats.gamma.cdf(x, shape, loc=loc, scale=scale)
        cdf = np.clip(cdf, 1e-6, 1 - 1e-6)    # keep the normal quantile finite
        out[sel] = stats.norm.ppf(cdf)
    return out


monthly = pd.read_csv(os.path.join(DATA, 'monthly_precip.csv'))
frames = []

for sid, g in monthly.groupby('station_id'):
    label = g['station'].iloc[0]
    g = g.sort_values(['year', 'month'])
    idx = pd.PeriodIndex(
        pd.to_datetime(dict(year=g.year, month=g.month, day=1)), freq='M')
    s = pd.Series(g['prcp_mm'].values, index=idx)

    # Reindex onto a gap-free monthly axis. Months dropped by the completeness
    # screen stay NaN, so any accumulation window containing one is also NaN
    # rather than being silently treated as a dry month.
    full = pd.period_range(idx.min(), idx.max(), freq='M')
    s = s.reindex(full)

    df = pd.DataFrame({'prcp_mm': s})
    df['station_id'] = sid
    df['station'] = label
    df['year'] = full.year
    df['month'] = full.month
    for k in SCALES:
        acc = s.rolling(k, min_periods=k).sum()
        df['spi%d' % k] = spi_series(acc, pd.Series(full.month, index=full))
    frames.append(df.reset_index(drop=True))

spi = pd.concat(frames, ignore_index=True)
spi['season'] = spi['month'].map(SEASONS)
spi.to_csv(os.path.join(DATA, 'spi.csv'), index=False)
print('SPI written for %d station-months' % spi['spi12'].notna().sum())

# ---- drought events: runs of consecutive months at or below the threshold ----
events = []
for (sid, label), g in spi.groupby(['station_id', 'station']):
    for k in SCALES:
        col = 'spi%d' % k
        v = g[col].values
        in_drought, start, acc = False, None, []
        for i, x in enumerate(v):
            wet = pd.isna(x) or x > DROUGHT_THRESHOLD
            if not wet and not in_drought:
                in_drought, start, acc = True, i, [x]
            elif not wet and in_drought:
                acc.append(x)
            elif wet and in_drought:
                events.append(dict(
                    station_id=sid, station=label, scale=k,
                    start_year=int(g['year'].iloc[start]), start_month=int(g['month'].iloc[start]),
                    end_year=int(g['year'].iloc[i - 1]), end_month=int(g['month'].iloc[i - 1]),
                    duration_months=len(acc),
                    severity=float(-np.sum(acc)),        # accumulated deficit
                    peak_spi=float(np.min(acc)),
                    intensity=float(-np.mean(acc))))
                in_drought = False
        if in_drought:
            events.append(dict(
                station_id=sid, station=label, scale=k,
                start_year=int(g['year'].iloc[start]), start_month=int(g['month'].iloc[start]),
                end_year=int(g['year'].iloc[len(v) - 1]), end_month=int(g['month'].iloc[len(v) - 1]),
                duration_months=len(acc), severity=float(-np.sum(acc)),
                peak_spi=float(np.min(acc)), intensity=float(-np.mean(acc))))

ev = pd.DataFrame(events)
ev['severity_class'] = ev['peak_spi'].map(classify)
ev.to_csv(os.path.join(DATA, 'drought_events.csv'), index=False)

# ---- summaries ----
by_station = []
for (sid, label), g in spi.groupby(['station_id', 'station']):
    for k in SCALES:
        col, e = 'spi%d' % k, ev[(ev.station_id == sid) & (ev.scale == k)]
        n = g[col].notna().sum()
        by_station.append(dict(
            station_id=sid, station=label, scale=k, months_with_spi=int(n),
            pct_months_in_drought=round(100 * (g[col] <= DROUGHT_THRESHOLD).sum() / n, 1) if n else np.nan,
            n_events=len(e),
            mean_duration_months=round(e.duration_months.mean(), 1) if len(e) else np.nan,
            max_duration_months=int(e.duration_months.max()) if len(e) else 0,
            mean_severity=round(e.severity.mean(), 2) if len(e) else np.nan,
            max_severity=round(e.severity.max(), 2) if len(e) else np.nan,
            min_spi=round(g[col].min(), 2) if n else np.nan))
bs = pd.DataFrame(by_station)
bs.to_csv(os.path.join(DATA, 'summary_by_station.csv'), index=False)

by_season = []
for (season, k), g in spi.melt(
        id_vars=['station_id', 'station', 'year', 'month', 'season'],
        value_vars=['spi%d' % k for k in SCALES],
        var_name='scale_col', value_name='spi').assign(
        scale=lambda d: d.scale_col.str.replace('spi', '').astype(int)
        ).groupby(['season', 'scale']):
    n = g['spi'].notna().sum()
    by_season.append(dict(
        season=season, scale=k, months_with_spi=int(n),
        pct_months_in_drought=round(100 * (g['spi'] <= DROUGHT_THRESHOLD).sum() / n, 1) if n else np.nan,
        mean_spi=round(g['spi'].mean(), 3),
        pct_moderate=round(100 * g['spi'].between(-1.5, -1.0, inclusive='right').sum() / n, 1) if n else np.nan,
        pct_severe=round(100 * g['spi'].between(-2.0, -1.5, inclusive='right').sum() / n, 1) if n else np.nan,
        pct_extreme=round(100 * (g['spi'] <= -2.0).sum() / n, 1) if n else np.nan))
pd.DataFrame(by_season).to_csv(os.path.join(DATA, 'summary_by_season.csv'), index=False)

print()
print('Drought events found (SPI <= %.1f):' % DROUGHT_THRESHOLD)
print(ev.groupby('scale').agg(events=('duration_months', 'size'),
                              mean_duration=('duration_months', 'mean'),
                              longest=('duration_months', 'max'),
                              worst_severity=('severity', 'max')).round(2).to_string())
print()
print('Percent of months in drought, SPI-12, by station:')
print(bs[bs.scale == 12][['station', 'pct_months_in_drought', 'n_events',
                          'max_duration_months', 'min_spi']].to_string(index=False))
