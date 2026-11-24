"""
@author: jbrondex

Description:
------------
Generic sensitivity analysis of the progressive-loading prognostic runs (FS solved, no GM flux on lateral BC,
parameterized friction) to the friction / cavitation parameters C, dc, h0 and As.

Any number of simulations can be compared: simply add/remove entries in the SIMULATIONS list below.
For each simulation, the output file name is rebuilt from its parameters following the convention
    MyBirch_C{C}_dc{dc}m_h0{h0}m_As{As}_Start{start}May_StopLoad{stopload}May_Stop{stop}May__bed.dat
e.g. C=0.29, dc=5, h0=0.3, As=37000 -> MyBirch_C029_dc5m_h003m_As37000_Start16May_StopLoad23May_Stop31May__bed.dat
A 'filename' key can be given to bypass this convention for files that do not follow it.

Figures produced:
    Fig 1 : mean basal velocity <ub> vs time (+ radar)
    Fig 2 : horizontal surface velocity at the 4 time-lapse points (one panel per point) vs observations
    Fig 3 : mean front velocity projected in radar line-of-sight vs radar
    Fig 4 : fraction of grounded area vs time
    Fig 5 : surface velocity at selected time-lapse point(s) + front velocity vs observations
"""

################################################################################
# import #####
################################################################################
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from scipy.interpolate import LinearNDInterpolator


################################################################################
# USER INPUTS #####
################################################################################
###Directory containing the simulation outputs
PATHROOT_SIMU = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/6_PRONO_PROGRESSIVELOAD_FS_NoGMLat_PARAMETERIZEDFRICTION/GlacierOut/')
PATHROOT_OBS = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/Obs/')

###Default values: used for any key not given in a simulation entry
DEFAULTS = dict(
    C=0.3,          # friction coefficient [-]
    dc=5.0,         # [m]
    h0=0.3,         # [m]
    As=37000,       # [m a-1 MPa-3]
    start_day=16,   # day of May 2025 at which the simu starts
    stopload_day=23,  # last day of loading (May 2025)
    stop_day=31,    # day of May 2025 at which the simu stops
    timestep=0.01 / 365,  # time step size IN YEARS
)

###List of simulations to compare. Only the keys that differ from DEFAULTS need to be given.
###Optional keys: 'label' (legend entry, otherwise built automatically from the parameters that vary),
###               'color', 'linewidth', 'filename' (overrides the file name convention), 'pathroot'.
SIMULATIONS = [
    dict(C=0.26),
    dict(C=0.28),
    dict(C=0.29),
    dict(C=0.30),
    # dict(C=0.25, dc=15, h0=0.1),
    # dict(C=0.25, dc=7,  h0=0.8),
    # dict(C=0.25, As=20000),
    # dict(label='RS As=12000', filename='MyBirch_RS_C03_dc5m_tsp001d_Start16May_StopLoad23May_Stop31May_As12000__bed.dat'),
]

###Optionally add the reference simulation with no friction over detached areas
INCLUDE_NOFRICTION_REF = True
PATHROOT_SIMU_NOFRICTION = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/')
NOFRICTION_REF = dict(
    label='No friction over detached areas ($C$ = 0.3)',
    filename='MyBirch_RS_C03_dc5m_tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat',
    pathroot=PATHROOT_SIMU_NOFRICTION,
    color='grey',
    reference=True,  ###excluded from automatic legend/title/color building
)

###Plot options
SIMU_DURATION = 15   ###x-axis extent in days from the (earliest) simulation start
LINEWIDTH = 3.0
COLORMAP = 'viridis'  ###used when no 'color' is given and there are more than 10 simulations (tab10 otherwise)
YLIMS = dict(ub=(0, 40), timelapse=(0, 200), front=(0, 40), grounded=(0, 1), fig5=(0, 100))
###Ids of time-lapse points (among 1, 2, 3, 4) displayed in Fig 5
POINTS_TIMELAPSE_IDS_FIG5 = [3]

###Fixed data
COLLAPSE_TIME = pd.Timestamp("2025-05-28 15:30")
RADAR_POSITION = np.array([2628689, 1141335, 1794])  ##xradar, yradar, zradar
LAST_RELIABLE_RADAR = pd.Timestamp("2025-05-25 06:00")  ###Radar data not reliable after 6 am on 25 May 25
POINTS_TIMELAPSE = np.array([[2630784, 1139171], [2630758, 1139194], [2630712, 1139233], [2630681, 1139269]])
COLORS_POINTS_TIMELAPSE = ['blue', 'orange', 'g', 'r', '#cc00cc']

COL_NAME = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'Normal Stress', 'theta_dagger', 'theta',
            'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
BC_FRONT, BC_BED, BC_TOP = 2, 3, 4

###Display names of parameters (for automatic legends/titles)
PARAM_DISPLAY = {
    'C': (r'$C$', ''),
    'dc': (r'$d_c$', ' m'),
    'h0': (r'$h_0$', ' m'),
    'As': (r'$A_s$', r' m a$^{-1}$ MPa$^{-3}$'),
    'stopload_day': ('end of loading', ' May'),
}


################################################################################
# FUNCTIONS #####
################################################################################
def fmt_noDot(x):
    """0.3 -> '03', 0.25 -> '025', 5 -> '5' (convention used for C and h0 in file names)."""
    s = f"{x:.6f}".rstrip('0').rstrip('.')
    return s.replace('.', '')


def fmt_dc(dc):
    """Convention used for dc in file names: 5 -> '5', 0.5 -> '05', 1.5 -> '1dot5'."""
    if float(dc).is_integer():
        return f"{int(dc)}"
    if dc < 1:
        return fmt_noDot(dc)
    x, y = str(dc).split('.')
    return f"{x}dot{y}"


def build_filename(p):
    return (f"MyBirch_C{fmt_noDot(p['C'])}_dc{fmt_dc(p['dc'])}m_h0{fmt_noDot(p['h0'])}m_As{int(p['As'])}"
            f"_Start{p['start_day']}May_StopLoad{p['stopload_day']}May_Stop{p['stop_day']}May__bed.dat")


def fmt_value(x):
    return f"{x:g}"


def complete_simulations(simulations):
    """Merge each simulation with DEFAULTS, build file name, label and color."""
    sims = [{**DEFAULTS, **s} for s in simulations]
    ###Reference simulations (e.g. full friction) are not used to build legends/title/colors
    param_sims = [s for s in sims if not s.get('reference', False)] or sims
    ###Parameters that differ between simulations -> shown in legend; the others -> shown in title
    varying = [k for k in PARAM_DISPLAY if len({s[k] for s in param_sims}) > 1]
    common = [k for k in PARAM_DISPLAY if k not in varying and k != 'stopload_day']
    if len(param_sims) <= 10:
        cmap_colors = plt.get_cmap('tab10').colors
    else:
        cmap_colors = plt.get_cmap(COLORMAP)(np.linspace(0, 0.95, len(param_sims)))
    i_color = 0
    for s in sims:
        s.setdefault('pathroot', PATHROOT_SIMU)
        s.setdefault('filename', build_filename(s))
        if 'color' not in s:
            s['color'] = cmap_colors[i_color % len(cmap_colors)]
            i_color += 1
        s.setdefault('linewidth', LINEWIDTH)
        if 'label' not in s:
            keys = varying if varying else ['C']
            s['label'] = ', '.join(f"{PARAM_DISPLAY[k][0]} = {fmt_value(s[k])}{PARAM_DISPLAY[k][1]}" for k in keys)
    title = ', '.join(f"{PARAM_DISPLAY[k][0]} = {fmt_value(sims[0][k])}{PARAM_DISPLAY[k][1]}" for k in common)
    return sims, title


def process_simulation(sim):
    """Load a bed .dat file and compute all time series needed for the figures."""
    filepath = Path(sim['pathroot']).joinpath(sim['filename'])
    print('Opening file:', filepath.name)
    df = pd.read_csv(filepath, names=COL_NAME, sep=r'\s+')
    t0 = pd.Timestamp(f"2025-05-{sim['start_day']} 00:00")
    df['Date'] = t0 + pd.to_timedelta(sim['timestep'] * df['Timestep'] * 365, unit='D')

    front = df[df['BC'] == BC_FRONT].copy()
    bed = df[df['BC'] == BC_BED].copy()
    top = df[df['BC'] == BC_TOP].copy()
    out = {}

    ####MEAN BASAL VELOCITY: tangential part of velocity (un should be zero but there is some residual)
    un = bed['u'] * bed['nx'] + bed['v'] * bed['ny'] + bed['w'] * bed['nz']
    bed['ub'] = np.sqrt((bed['u'] - un * bed['nx']) ** 2 + (bed['v'] - un * bed['ny']) ** 2 + (bed['w'] - un * bed['nz']) ** 2)
    out['mean_ub'] = bed.groupby('Date')['ub'].mean() / 365  ### m/day

    ####GROUNDED AREA FRACTION
    total_area = bed.groupby('Date')['nodearea'].sum()
    grounded_area = bed[bed['GM'] == 1].groupby('Date')['nodearea'].sum()
    out['grounded_ratio'] = (grounded_area / total_area).fillna(0.0)

    ####FRONT VELOCITY PROJECTED ON RADAR LINE OF SIGHT
    d = front[['X', 'Y', 'Z']].values - RADAR_POSITION
    r_norm = np.linalg.norm(d, axis=1)
    front['uf_radarproj'] = np.einsum('ij,ij->i', front[['u', 'v', 'w']].values, d) / r_norm
    out['mean_uf_radar'] = front.groupby('Date')['uf_radarproj'].mean().abs() / 365  ### m/day

    ####HORIZONTAL SURFACE VELOCITY INTERPOLATED AT TIME-LAPSE POINTS
    rows = []
    for t, group in top.groupby('Date'):
        interp = LinearNDInterpolator(group[['X', 'Y']].values, group[['u', 'v']].values)
        uv = interp(POINTS_TIMELAPSE)
        for pid, (u_i, v_i) in enumerate(uv, start=1):
            rows.append({'Date': t, 'Point_id': pid, 'u_interp_norm': np.hypot(u_i, v_i) / 365})  ### m/day
    out['timelapse'] = pd.DataFrame(rows)
    return out


def load_observations():
    ######## Radar Data (BEWARE: in mm/day) #########
    col_radar = ['Time', 'Fullfront_6h', 'Middlefront_6h', 'Upperfront_6h', 'Fullfront_24h', 'Middlefront_24h',
                 'Upperfront_24h', 'Fullfront_RWB', 'Middlefront_RWB', 'Upperfront_RWB', 'FullMean']
    radar = pd.read_csv(PATHROOT_OBS.joinpath('glacier_front_RadarData_mmperday.csv'), header=0, names=col_radar)
    radar['Date'] = pd.to_datetime(radar['Time'], format="%d.%m.%y %H:%M")
    radar = radar[radar['Date'] < LAST_RELIABLE_RADAR].copy()
    ###Mean of the three processing methods (6h vs 24h vs RWB) for the full front, converted to m/day
    radar['mean_Fullfront'] = radar[['Fullfront_6h', 'Fullfront_24h', 'Fullfront_RWB']].mean(axis=1, skipna=True).abs() / 1000

    ####### TimeLapse Velocities ####
    col_tl = ['X0', 'Y0', 't0', 't1', 'horizontal_velo', 'sigma']
    tl = pd.read_csv(PATHROOT_OBS.joinpath('VelocitiesFromTimelapse.csv'), header=0, names=col_tl)
    tl['t0'] = pd.to_datetime(tl['t0'], utc=True)
    tl['t1'] = pd.to_datetime(tl['t1'], utc=True)
    tl['t_mid'] = tl['t0'] + (tl['t1'] - tl['t0']) / 2
    ###Model dates are tz-naive: drop tz info so that both can share the same axis
    tl['t_mid'] = tl['t_mid'].dt.tz_localize(None)
    return radar, tl


def timelapse_obs_at_point(tl, pid):
    X, Y = POINTS_TIMELAPSE[pid - 1]
    return tl[(tl['X0'] == X) & (tl['Y0'] == Y)]


def format_time_axis(ax, t0, tf, ylim=None):
    ax.set_xlim(t0, tf)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
    ax.tick_params(labelsize=20)
    ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
    ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
    for lab in ax.get_xticklabels():
        lab.set_rotation(30)
        lab.set_ha('right')


def add_collapse_line(ax, fontsize=22):
    ax.axvline(COLLAPSE_TIME, linestyle='--', linewidth=2.7, color='r')
    ax.text(COLLAPSE_TIME, 0.98, 'Collapse', fontsize=fontsize, color='r', rotation=90, va='top', ha='right',
            transform=ax.get_xaxis_transform())


def new_fig(ylabel, title, figsize=(15, 12)):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.set_ylabel(ylabel, fontsize=26)
    ax.set_title(title, fontsize=24, weight='bold')
    return fig, ax


################################################################################
# MAIN #####
################################################################################
if __name__ == "__main__":
    for pid in POINTS_TIMELAPSE_IDS_FIG5:
        if pid not in {1, 2, 3, 4}:
            raise ValueError(f"Invalid time-lapse point id {pid}: should be 1, 2, 3 or 4.")

    all_simulations = SIMULATIONS + ([NOFRICTION_REF] if INCLUDE_NOFRICTION_REF else [])
    sims, common_title = complete_simulations(all_simulations)
    Df_Radar, Df_Timelapse = load_observations()

    ###Process all simulations (skip missing files with a warning)
    results = []
    for sim in sims:
        try:
            results.append((sim, process_simulation(sim)))
        except FileNotFoundError:
            print(f"WARNING: file not found, simulation skipped: {sim['filename']}")
    if not results:
        raise RuntimeError("None of the requested simulations could be loaded.")

    t0 = min(pd.Timestamp(f"2025-05-{s['start_day']} 00:00") for s, _ in results)
    tf = t0 + pd.Timedelta(days=SIMU_DURATION)
    suptitle = f"Common parameters: {common_title}" if common_title else None

    ###~~~~~~~ FIG 1 : mean ub ~~~~~~~
    fig1, ax1 = new_fig(r'$\langle u_\mathrm{b} \rangle$ [m d$^{-1}$]', 'Mean basal velocity')
    ###~~~~~~~ FIG 2 : one panel per time-lapse point ~~~~~~~
    fig2, axs2 = plt.subplots(2, 2, figsize=(20, 14), constrained_layout=True, sharex=True, sharey=True)
    axs2 = axs2.ravel()
    ###~~~~~~~ FIG 3 : front velocity in radar direction ~~~~~~~
    fig3, ax3 = new_fig(r'$\langle u_\mathrm{f} \rangle$ [m d$^{-1}$]', 'Velocity at front vs radar')
    ###~~~~~~~ FIG 4 : grounded area fraction ~~~~~~~
    fig4, ax4 = new_fig('Grounded Area', 'Fraction of grounded area')
    ###~~~~~~~ FIG 5 : surface velocity at selected point(s) + front velocity ~~~~~~~
    fig5, ax5 = new_fig(r'Ice velocity [m d$^{-1}$]', common_title)

    for sim, res in results:
        kw = dict(color=sim['color'], linewidth=sim['linewidth'])
        ax1.plot(res['mean_ub'].index, res['mean_ub'], label=sim['label'], linestyle='-', **kw)
        ax3.plot(res['mean_uf_radar'].index, res['mean_uf_radar'], label=sim['label'], linestyle='-', **kw)
        ax4.plot(res['grounded_ratio'].index, res['grounded_ratio'], label=sim['label'], linestyle='-', **kw)
        for pid, ax in enumerate(axs2, start=1):
            sub = res['timelapse'][res['timelapse']['Point_id'] == pid]
            ax.plot(sub['Date'], sub['u_interp_norm'], label=sim['label'], linestyle='-', **kw)
        for pid in POINTS_TIMELAPSE_IDS_FIG5:
            sub = res['timelapse'][res['timelapse']['Point_id'] == pid]
            ax5.plot(sub['Date'], sub['u_interp_norm'], linestyle='-', **kw)
        ax5.plot(res['mean_uf_radar'].index, res['mean_uf_radar'], linestyle='--',
                 color=sim['color'], linewidth=sim['linewidth'] * 0.8)

    ###Observations
    radar_kw = dict(color='k', linestyle='--', linewidth=LINEWIDTH * 1.6, label='Radar')
    for ax in (ax1, ax3):
        ax.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'], **radar_kw)
    for pid, ax in enumerate(axs2, start=1):
        obs = timelapse_obs_at_point(Df_Timelapse, pid)
        ax.scatter(obs['t_mid'], obs['horizontal_velo'], color='k', s=25, label='Time-lapse obs.', zorder=5)
        ax.set_title(f'Point {pid}', fontsize=22, weight='bold')
        ax.set_ylabel(r'$|u_\mathrm{s,h}|$ [m d$^{-1}$]', fontsize=20)
    ax5.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'], color='k', linestyle='--', linewidth=LINEWIDTH * 2)
    marker_colors = COLORS_POINTS_TIMELAPSE[-1:] + COLORS_POINTS_TIMELAPSE
    for pid, col in zip(POINTS_TIMELAPSE_IDS_FIG5, marker_colors):
        obs = timelapse_obs_at_point(Df_Timelapse, pid)
        ax5.scatter(obs['t_mid'], obs['horizontal_velo'], color=col, s=20, zorder=5)

    ###Axes formatting + collapse line
    format_time_axis(ax1, t0, tf, YLIMS['ub'])
    format_time_axis(ax3, t0, tf, YLIMS['front'])
    format_time_axis(ax4, t0, tf, YLIMS['grounded'])
    format_time_axis(ax5, t0, tf, YLIMS['fig5'])
    for ax in axs2:
        format_time_axis(ax, t0, tf, YLIMS['timelapse'])
    for ax in (ax1, ax3, ax4, ax5, *axs2):
        add_collapse_line(ax)

    ###Legends
    for ax, loc in ((ax1, 'upper left'), (ax3, 'upper left'), (ax4, 'lower left')):
        ax.legend(fontsize=18, loc=loc)
    axs2[0].legend(fontsize=14, loc='upper left')
    if suptitle:
        for fig in (fig1, fig2, fig3, fig4):
            fig.suptitle(suptitle, fontsize=20)

    ###Fig 5: style legend (what) + color legend (which simu)
    style_handles = [
        Line2D([], [], color='grey', linestyle='-', linewidth=LINEWIDTH, label='Modelled velocity (surface point)'),
        Line2D([], [], color=marker_colors[0], marker='o', linestyle='none', label='Time-lapse observations (surface point)'),
        Line2D([], [], color='k', linestyle='--', linewidth=LINEWIDTH * 2, label='Radar observations (front)'),
        Line2D([], [], color='grey', linestyle='--', linewidth=LINEWIDTH * 0.8, label='Modelled front velocity'),
    ]
    sim_handles = [Line2D([], [], color=s['color'], linewidth=s['linewidth'], label=s['label']) for s, _ in results]
    leg5a = ax5.legend(handles=style_handles, fontsize=18, loc='upper left')
    ax5.add_artist(leg5a)
    ax5.legend(handles=sim_handles, fontsize=18, loc='center left')

    plt.show()
