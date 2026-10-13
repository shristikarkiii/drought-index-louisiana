"""
Figures for the Louisiana drought analysis.

    python make_figures.py      ->  figures/*.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, FIGS = os.path.join(HERE, 'data'), os.path.join(HERE, 'figures')
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.edgecolor': '#6b6b6b', 'axes.linewidth': .8,
                     'axes.grid': True, 'grid.color': '#e6e3dd', 'grid.linewidth': .7,
                     'axes.axisbelow': True, 'figure.dpi': 200})

DRY, WET = '#8c3b1e', '#1f5f8b'
CMAP = LinearSegmentedColormap.from_list(
    'spi', ['#6b2d12', '#b5561f', '#e2a76b', '#f5efe6', '#9dc3dd', '#3d7fae', '#17465f'])

spi = pd.read_csv(os.path.join(DATA, 'spi.csv'))
ev = pd.read_csv(os.path.join(DATA, 'drought_events.csv'))
bs = pd.read_csv(os.path.join(DATA, 'summary_by_station.csv'))
bseas = pd.read_csv(os.path.join(DATA, 'summary_by_season.csv'))

spi['t'] = spi['year'] + (spi['month'] - 1) / 12.0
# north to south, so the heatmap reads geographically
ORDER = ['Plain Dealing (NW)', 'Homer (N)', 'Winnsboro (NE)', 'Alexandria (Central)',
         'Lafayette Regional AP (SW)', 'New Iberia (S Central)',
         'Covington (SE)', 'New Orleans Audubon (SE Coastal)']


def save(fig, name):
    fig.savefig(os.path.join(FIGS, name), bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('  wrote', name)


# 1 -- the headline: SPI-12 for every station through time
piv = spi.pivot_table(index='station', columns='t', values='spi12')
piv = piv.reindex([s for s in ORDER if s in piv.index])
fig, ax = plt.subplots(figsize=(11, 3.6))
im = ax.pcolormesh(piv.columns.values, np.arange(len(piv)), piv.values,
                   cmap=CMAP, norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3),
                   shading='nearest')
ax.set_yticks(np.arange(len(piv)))
ax.set_yticklabels(piv.index, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel('Year')
ax.grid(False)
ax.set_title('12-month SPI across the Louisiana station network, 1950–2026',
             fontsize=11.5, fontweight='bold', loc='left', pad=22)
ax.text(0, 1.015,
        'Brown is drier than the local long-term normal for that month, blue is wetter. '
        'White marks months the completeness screen rejected.',
        transform=ax.transAxes, fontsize=8.4, color='#5d6360', va='bottom')
cb = fig.colorbar(im, ax=ax, fraction=.025, pad=.012)
cb.set_label('SPI-12 (standard deviations)', fontsize=8.5)
save(fig, 'fig1_spi12_heatmap.png')

# 2 -- drought frequency, duration and severity by station
sub = bs[bs.scale == 12].set_index('station').reindex([s for s in ORDER if s in bs.station.values])
fig, axes = plt.subplots(1, 3, figsize=(11, 3.3))
for ax, col, lab, c in [
        (axes[0], 'pct_months_in_drought', 'Months in drought (%)', DRY),
        (axes[1], 'mean_duration_months', 'Mean event duration (months)', '#a8632c'),
        (axes[2], 'max_severity', 'Worst event severity (SPI-months)', '#6b2d12')]:
    ax.barh(sub.index[::-1], sub[col][::-1], color=c, height=.68)
    ax.set_xlabel(lab, fontsize=8.6)
    ax.tick_params(labelsize=7.6)
    if ax is not axes[0]:
        ax.set_yticklabels([])
axes[0].axvline(15.87, color='#333', lw=1, ls='--')
axes[0].text(15.87, -.9, ' expected 15.9%', fontsize=7, color='#333', va='top')
fig.suptitle('Drought characteristics by station, 12-month SPI',
             fontsize=11.5, fontweight='bold', x=.09, ha='left', y=1.03)
save(fig, 'fig2_by_station.png')

# 3 -- season
fig, ax = plt.subplots(figsize=(6.6, 3.4))
order = ['Winter (DJF)', 'Spring (MAM)', 'Summer (JJA)', 'Autumn (SON)']
w, scales = .26, [3, 6, 12]
for i, k in enumerate(scales):
    d = bseas[bseas.scale == k].set_index('season').reindex(order)
    ax.bar(np.arange(4) + (i - 1) * w, d['pct_months_in_drought'], width=w,
           label='SPI-%d' % k, color=['#e2a76b', '#b5561f', '#6b2d12'][i])
ax.axhline(15.87, color='#333', lw=1, ls='--')
ax.text(3.45, 16.4, 'expected 15.9%', fontsize=7.4, color='#333', ha='right')
ax.set_xticks(np.arange(4)); ax.set_xticklabels(order, fontsize=8.4)
ax.set_ylabel('Months meeting drought threshold (%)')
ax.legend(frameon=False, fontsize=8.4)
ax.set_title('Drought incidence by season and accumulation period',
             fontsize=11.5, fontweight='bold', loc='left', pad=11)
save(fig, 'fig3_by_season.png')

# 4 -- statewide mean, to show the shared episodes
fig, ax = plt.subplots(figsize=(11, 3.1))
st = spi.groupby('t')['spi12'].mean()
ax.fill_between(st.index, 0, st.where(st < 0), color=DRY, lw=0)
ax.fill_between(st.index, 0, st.where(st > 0), color=WET, lw=0)
ax.axhline(-1, color='#333', lw=.9, ls='--')
ax.text(st.index.min() + .4, -1.06, 'drought threshold, SPI = −1',
        fontsize=7.6, color='#333', va='top')
ax.set_ylabel('Mean SPI-12')
ax.set_xlabel('Year')
ax.set_title('Statewide mean 12-month SPI: shared drought episodes',
             fontsize=11.5, fontweight='bold', loc='left', pad=11)
save(fig, 'fig4_statewide.png')

# 5 -- distribution check, evidence the fit behaves
fig, ax = plt.subplots(figsize=(5.2, 3.3))
v = spi['spi12'].dropna()
ax.hist(v, bins=48, density=True, color='#c7bfb2', edgecolor='white', linewidth=.5)
x = np.linspace(-4, 4, 200)
ax.plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color='#6b2d12', lw=1.6,
        label='standard normal')
ax.axvline(-1, color='#333', lw=1, ls='--')
ax.set_xlabel('SPI-12'); ax.set_ylabel('Density')
ax.legend(frameon=False, fontsize=8.4)
ax.set_title('SPI-12 distribution against the standard normal',
             fontsize=10.5, fontweight='bold', loc='left', pad=10)
ax.text(0, -.30, 'mean %.3f, sd %.3f, n = %d' % (v.mean(), v.std(), len(v)),
        transform=ax.transAxes, fontsize=7.8, color='#5d6360')
save(fig, 'fig5_distribution.png')

print()
worst = ev[ev.scale == 12].nlargest(6, 'severity')[
    ['station', 'start_year', 'start_month', 'duration_months', 'severity', 'peak_spi']]
print('Six most severe 12-month-SPI droughts on record:')
print(worst.to_string(index=False))
