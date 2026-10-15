"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of velocities (at front and at selected points) and grounded area over time to compare results obtained for
various values of C and dc using the most recent implementation of effective pressure in the rate-and-state friction law (July 2026)
where effective pressure is computed from the true normal stress at bed

"""
from cProfile import label

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
# Create legend handles manually, for example:
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path

import numpy as np
import re
import os
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
from matplotlib import ticker
from scipy.interpolate import griddata
from matplotlib.ticker import MultipleLocator
import matplotlib.dates as mdates

# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-3, 3))

###Defined Colors optimized for color-blind people:
RedStar = [162 / 255, 20 / 255, 47 / 255]
RedBackGround = [233 / 255, 201 / 255, 207 / 255]
LightBlueGreyBackGround = [194 / 255, 209 / 255, 219 / 255]
BlueGreyBackGround = [147/ 255, 174 / 255, 191 / 255]
LightBrownBackGround = [223/ 255, 211 / 255, 205 / 255]
VeryLightBrown = [255 / 255, 199 / 255, 127 / 255]
LightBrown = [229 / 255, 143 / 255, 91 / 255]
Brown = [140 / 255, 88 / 255, 56 / 255]
DarkBrown = [64 / 255, 40 / 255, 25 / 255]
VeryDarkBrown = [13 / 255, 8 / 255, 5 / 255]

####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN PART OF THE CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####


###############################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    #### SCRIPT PARAMETERS ###
    ###Parameter of Simu
    Cvalues= [0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39]
    ###Corresponding colors in this order: black, blue, orange, green, brown, purple
    colors = ['k','#377eb8', '#ff7f00', '#4daf4a', '#a65628', '#984ea3','#e41a1c',    # red
    '#f781bf',    # pink
    '#999999',    # gray
    '#dede00'     # yellow (olive-like, better visibility than pure yellow)
]
    LineSize = 2.4

    Year = '2025'
    StartSimuDay = 16  ####Day of May 2025 corresponding to simulation start
    SimuDuration = 15  ###Simu duration in days
    EndOfLoadingDays = 23
    t0 = pd.Timestamp("2025-05-{} 00:00".format(StartSimuDay))
    tf = t0 + pd.Timedelta(days=SimuDuration)
    CollapseTime = pd.Timestamp("2025-05-28 15:30")

    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'Normal Stress', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    Pathroot_July26 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')

    ###TREAT OBSERVATIONS###
    ######## Radar Data #########
    RadarPosition = [2628689, 1141335,1794] ##xradar, yradar, zradar
    LastReliableRadarData = pd.Timestamp("2025-05-25 06:00") ###Radar data is not reliable anymore after 6 am on the 25th of May 25
    Col_Name_Radar = ['Time','Fullfront_6h','Middlefront_6h','Upperfront_6h','Fullfront_24h','Middlefront_24h','Upperfront_24h','Fullfront_RWB','Middlefront_RWB','Upperfront_RWB','FullMean']
    Pathroot_Radar = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/Obs/.')
    Filename_Radar = 'glacier_front_RadarData_mmperday.csv'
    ###BE WARE : Radar data are in mm/day
    Df_Radar = pd.read_csv(Pathroot_Radar.joinpath(Filename_Radar), header=0, names=Col_Name_Radar)
    ###Convert time of csv file in panda date
    Df_Radar['Date'] = pd.to_datetime(Df_Radar["Time"], format="%d.%m.%y %H:%M")
    ###Cut all data that are not reliable (i.e. aquired later than last reliable data)
    Df_Radar =  Df_Radar[Df_Radar['Date']<LastReliableRadarData]
    ###Calculate a mean between the three processing method (6h vs 24h vs RWB) for the full front
    Df_Radar["mean_Fullfront"] = Df_Radar[["Fullfront_6h", "Fullfront_24h", "Fullfront_RWB"]].mean(axis=1, skipna=True)

    ####### TimeLapse Velocities ####
    ###Coordinates of surface points at which horizontal velocity is inferred
    Point1_Timelapse = [2630784,1139171]
    Point2_Timelapse = [2630758,1139194]
    Point3_Timelapse = [2630712,1139233]
    Point4_Timelapse = [2630681,1139269]
    Colors_Points_Timelapse = ['blue', 'orange','g','r']
    Col_Name_Timelapse = ['X0','Y0','t0','t1','horizontal_velo','sigma']
    Pathroot_Timelapse = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/Obs/.')
    Filename_Timelapse = 'VelocitiesFromTimelapse.csv'
    Df_Timelapse = pd.read_csv(Pathroot_Timelapse.joinpath(Filename_Timelapse), header=0, names=Col_Name_Timelapse)
    ###Get the mean time between two images used to compute velocities
    Df_Timelapse['t0'] = pd.to_datetime(Df_Timelapse['t0'], utc=True)
    Df_Timelapse['t1'] = pd.to_datetime(Df_Timelapse['t1'], utc=True)
    Df_Timelapse['t_mid'] = Df_Timelapse['t0'] + (Df_Timelapse['t1'] - Df_Timelapse['t0']) / 2

    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=20)
    ax.set_title('Mean basal velocity [m/day]', fontsize=21, weight='bold')
    ax.tick_params(labelsize=24)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)

    ###~~~~~~~FIG 2 : Horizontal velocity at selected surface points + mean ub ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig2, ax2 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax2.set_ylabel(r'velocity [m/day]', fontsize=30)
    ax2.set_yscale("linear")
    ax2.tick_params(labelsize=24)  # fontsize of the tick labels
    ax2.grid(True)
    ax2.grid(alpha=0.5)
    ax2.set_title('Horizontal velocities versus timelapse velocities', fontsize=30, weight='bold')

    ###~~~~~~~FIG 3 : Front velocity in interferometric radar direction ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig3, ax3 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax3.set_ylabel(r'$\langle u_\mathrm{f} \rangle$ [m/day]', fontsize=30)
    ax3.set_yscale("linear")
    ax3.tick_params(labelsize=24)  # fontsize of the tick labels
    ax3.grid(True)
    ax3.grid(alpha=0.5)
    ax3.set_title('Velocity at front vs radar', fontsize=30, weight='bold')

    ###~~~~~~~FIG 4 : Grounded area fraction ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig4, ax4 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax4.set_ylabel(r'Grounded Area', fontsize=30)
    ax4.set_yscale("linear")
    ax4.tick_params(labelsize=24)  # fontsize of the tick labels
    ax4.grid(True)
    ax4.grid(alpha=0.5)
    ax4.set_title('Fraction of grounded area', fontsize=30, weight='bold')

    ### START OF SCRIPT###
    ###Loop on each C considered
    for i, (Cvalue, blue_col) in enumerate(zip(Cvalues, colors)):
        print('Dealing with C value :', Cvalue)
        ####Rearrangement of parameters to get corresponding file names
        s = f"{Cvalue:.2f}"
        if len(s.rstrip('0').split('.')[1]) == 1:
            str_Cvalue = f"{int(Cvalue * 10):02d}"
        elif len(s.rstrip('0').split('.')[1]) == 2:
            str_Cvalue = f"{int(Cvalue * 100):03d}"

        if Cvalue == 0.3:
            dc_values = [2, 3, 4, 5]
        elif 0.31 <= Cvalue <= 0.34:
            dc_values =  [3]
        elif Cvalue == 0.35:
            dc_values = [2, 3]
        else:
            dc_values = [2]

        for dcvalue in dc_values:
            if dcvalue < 1:
                str_dcvalue = f"dc{int(dcvalue * 10):02d}m_"
            elif dcvalue.is_integer():
                str_dcvalue = f"dc{int(dcvalue):01d}m_"
            else:
                x, y = str(dcvalue).split(".")
                str_dcvalue = f"dc{x}dot{y}m_"

            Filename_July26 = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad{}May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue,EndOfLoadingDays)
            ###Time step size for prono simu with overload:
            Timestep_size = 0.01 / 365  ### In years !!
            ###load file corresponding to first simu as dataframe
            print('Opening file:', Filename_July26)
            Df_wholebed_July26 = pd.read_csv(Pathroot_July26.joinpath(Filename_July26), names=Col_Name, sep=r'\s+')
            ###########################################################################
            ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED END OF OVERLOAD DAY
            ###Create a column with time (in years !)
            Df_wholebed_July26['Time'] = Timestep_size * Df_wholebed_July26['Timestep']
            ###Convert absolute time in terms of date
            Df_wholebed_July26['Date'] = t0 + pd.to_timedelta(Df_wholebed_July26['Time'] * 365, unit='D')
            ###Create a dedicated dataframe for each boundary of interest
            Df_wholebed_front = Df_wholebed_July26.loc[Df_wholebed_July26['BC'] == 2].copy()
            Df_wholebed_bed = Df_wholebed_July26.loc[Df_wholebed_July26['BC'] == 3].copy()
            Df_wholebed_top = Df_wholebed_July26.loc[Df_wholebed_July26['BC'] == 4].copy()

            ####MEAN BASAL VELOCITY ####
            ###Calculate ub from velocity and normal vector
            ###Start by calculating un (should be zero but there is some residual)
            Df_wholebed_bed['un'] = Df_wholebed_bed['u'] * Df_wholebed_bed['nx'] + Df_wholebed_bed['v'] * Df_wholebed_bed['ny'] + Df_wholebed_bed['w'] * Df_wholebed_bed['nz']
            ###Deduce ub from u and un
            Df_wholebed_bed['ub'] = np.sqrt((Df_wholebed_bed['u'] - Df_wholebed_bed['un'] * Df_wholebed_bed['nx']) ** 2 + (Df_wholebed_bed['v'] - Df_wholebed_bed['un'] * Df_wholebed_bed['ny']) ** 2 + (Df_wholebed_bed['w'] - Df_wholebed_bed['un'] * Df_wholebed_bed['nz']) ** 2)
            ###Now calculate the mean basal velocity
            mean_ub_by_time = Df_wholebed_bed.groupby('Date')['ub'].mean()
            ###GROUNDED AREA FRACTION
            ### Calculate total area over the whole glacier for each time step
            Total_Area = Df_wholebed_bed.groupby('Date')['nodearea'].sum()
            ### Grounded area
            Total_GroundedArea = (Df_wholebed_bed[Df_wholebed_bed['GM'] == 1].groupby('Date')['nodearea'].sum())
            ### Fraction of grounded area over total area
            Ratio = Total_GroundedArea / Total_Area

            ####FRONT VELOCITY IN RADAR DIRECTION ####
            ###for each node of front, we calculate the vector from radar to node:
            dx = Df_wholebed_front['X'] - RadarPosition[0]
            dy = Df_wholebed_front['Y'] - RadarPosition[1]
            dz = Df_wholebed_front['Z'] - RadarPosition[2]
            ###for each node of front, we calculate norm of this vector:
            r_norm = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            ###Now do the projection through scalar product with unit vector
            Df_wholebed_front['uf_radarproj'] = (Df_wholebed_front['u'] * dx + Df_wholebed_front['v'] * dy + Df_wholebed_front['w'] * dz) / r_norm
            ###Now calculate the mean of all nodes at front for each Date
            mean_uf_radarproj_by_time = Df_wholebed_front.groupby('Date')['uf_radarproj'].mean()

            ####SURFACE VELOCITY####
            ###first calculate mean surface velocity
            Df_wholebed_top['us_norm'] = np.sqrt(Df_wholebed_top['u'] ** 2 + Df_wholebed_top['v'] ** 2 + Df_wholebed_top['w'] ** 2)
            ###Now calculate the mean norm of surface velocity for each Date
            mean_us_norm_by_time = Df_wholebed_top.groupby('Date')['us_norm'].mean()
            ###Now I interpolate the modeled horizontal velocity at the four points of the timelapse camera for each date
            results = []
            ##first I will do that for each date available
            for t, group in Df_wholebed_top.groupby('Date'):
                ###get x and y of all nodes
                nodes = group[['X', 'Y']].values
                # get horizontal velocities at nodes
                u_vals = group['u'].values
                v_vals = group['v'].values
                # do the interpolation at each of the four points of the timelapse camera
                # for id, Point in enumerate([Point1_Timelapse, Point2_Timelapse, Point3_Timelapse, Point4_Timelapse]):
                for id, Point in enumerate([Point3_Timelapse]):###just for one of the points
                    u_interp = griddata(nodes, u_vals, Point, method='linear')[0]
                    v_interp = griddata(nodes, v_vals, Point, method='linear')[0]
                    results.append({
                        'Date': t,
                        'Point_id': id + 1,  ##id starts at 1
                        'Xpt': Point[0],
                        'Ypt': Point[1],
                        'u_interp': u_interp,
                        'v_interp': v_interp
                    })
            Df_ModelInterpAtTimeLapsePoints = pd.DataFrame(results)
            ###Calculate norm of horizontal velocity
            Df_ModelInterpAtTimeLapsePoints['u_interp_norm'] = np.sqrt(Df_ModelInterpAtTimeLapsePoints['u_interp'] ** 2 + Df_ModelInterpAtTimeLapsePoints['v_interp'] ** 2)

            #################################################
            #####  FROM THERE WE CAN START THE FIGURES  #####
            #################################################
            ###Fig 1: Mean ub ~~~~~~~~~~~~~~~
            if dcvalue == 2:
                ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, label='C = {}'.format(Cvalue), color=blue_col, linestyle='-',linewidth=LineSize)
                if i == 0:##plot just once
                    ax.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs() / 1000, color='b', label='Radar', linestyle='--', linewidth=LineSize * 2)
                    ###Add vertical line for collapse day and annotate
                    ax.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
                    ax.text(CollapseTime, 40, 'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
            elif dcvalue == 3:
                if 0.31 <= Cvalue <= 0.34:
                    ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, label='C = {}'.format(Cvalue), color=blue_col, linestyle='--', linewidth=LineSize)
                else:
                    ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, color=blue_col,linestyle='--', linewidth=LineSize)
            elif dcvalue == 4:
                ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, color=blue_col,linestyle=':', linewidth=LineSize)
            elif dcvalue == 5:
                ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, color=blue_col,linestyle='-.', linewidth=LineSize)
            ###Dummy plots for legend regarding lateral B
            ax.set_ylim(0.0, 200.0)
            ax.set_xlim(t0, tf)
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig1.autofmt_xdate()

            ###Fig 2: Horizontal velocity at selected surface points + mean ub ~~~~~~~~~~~~~~~
            ### TimeLapse Velocities Modeled versus Observed
            # for l, Point in enumerate([Point1_Timelapse, Point2_Timelapse, Point3_Timelapse, Point4_Timelapse]):
            for l, Point in enumerate([Point3_Timelapse]):
                subset = Df_ModelInterpAtTimeLapsePoints[Df_ModelInterpAtTimeLapsePoints['Point_id'] == l + 1]
                subset_Obs = Df_Timelapse[(Df_Timelapse['X0'] == Point[0]) & (Df_Timelapse['Y0'] == Point[1])]
                ax2.scatter(subset_Obs['t_mid'], subset_Obs['horizontal_velo'], color=Colors_Points_Timelapse[l],s=30, label='Obs')
                if dcvalue == 2:
                    ax2.plot(subset['Date'], subset['u_interp_norm'] / 365, label='C = {}'.format(Cvalue),color=blue_col, linestyle='-', linewidth=LineSize)
                    if i == 0:
                        ### Radar Data
                        ax2.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs() / 1000, color='b', label='Radar',linestyle='--', linewidth=LineSize * 2)
                        ###Add vertical line for collapse day and annotate
                        ax2.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
                        ax2.text(CollapseTime + pd.Timedelta(hours=9, minutes=30), 190, 'Collapse', fontsize=26, color='r',rotation=90, va='top', ha='right')
                elif dcvalue == 3:
                    if 0.31 <= Cvalue <= 0.34:
                        ax2.plot(subset['Date'], subset['u_interp_norm'] / 365, label='C = {}'.format(Cvalue), color=blue_col, linestyle='--', linewidth=LineSize)
                    else:
                        ax2.plot(subset['Date'], subset['u_interp_norm'] / 365, color=blue_col, linestyle='--', linewidth=LineSize)
                elif dcvalue == 4:
                    ax2.plot(subset['Date'], subset['u_interp_norm'] / 365, color=blue_col, linestyle=':', linewidth=LineSize)
                elif dcvalue == 5:
                    ax2.plot(subset['Date'], subset['u_interp_norm'] / 365, color=blue_col, linestyle='-.', linewidth=LineSize)
            ###Dummy plots for legend regarding lateral B
            ax2.set_ylim(0.0, 200.0)
            ax2.set_xlim(t0, tf)
            ax2.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax2.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax2.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig2.autofmt_xdate()

            ###Fig 3: Velocity at front in radar direction~~~~~~~~~~~~~~~
            if dcvalue == 2:
                ax3.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365,label='C = {}'.format(Cvalue), color=blue_col, linestyle='-', linewidth=LineSize)
                if i == 0:
                    ax3.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs() / 1000, color='b', label='Radar', linestyle='--', linewidth=LineSize * 2)
                    ###Add vertical line for collapse day and annotate
                    ax3.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
                    ax3.text(CollapseTime, 40, 'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
            elif dcvalue == 3:
                if 0.31 <= Cvalue <= 0.34:
                    ax3.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365, color=blue_col, label='C = {}'.format(Cvalue),linestyle='--', linewidth=LineSize)
                else:
                    ax3.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365, color=blue_col, linestyle='--', linewidth=LineSize)
            elif dcvalue == 4:
                ax3.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365, color=blue_col, linestyle=':', linewidth=LineSize)
            elif dcvalue == 5:
                ax3.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365, color=blue_col, linestyle='-.', linewidth=LineSize)
            ###Dummy plots for legend regarding lateral B
            ax3.set_ylim(0.0, 40.0)
            ax3.set_xlim(t0, tf)
            ax3.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax3.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax3.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax3.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax3.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig3.autofmt_xdate()

            ###Fig 4: Grounded Area Fraction
            ### Plot ratio
            if dcvalue == 2:
                ax4.plot(Total_Area.index, Ratio, color=blue_col, label='C = {}'.format(Cvalue),linestyle='-', linewidth=LineSize)
                if i==0:
                ###Add vertical line for collapse day
                    ax4.axvline(CollapseTime,linestyle='--',linewidth=2.7, color= 'r')
                    ax4.text(CollapseTime, 40, 'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
            elif dcvalue == 3:
                if 0.31 <= Cvalue <= 0.34:
                    ax4.plot(Total_Area.index, Ratio, color=blue_col, label='C = {}'.format(Cvalue), linestyle='--', linewidth=LineSize)
                else:
                    ax4.plot(Total_Area.index, Ratio, color=blue_col, linestyle='--', linewidth=LineSize)
            elif dcvalue == 4:
                ax4.plot(Total_Area.index, Ratio, color=blue_col, linestyle=':', linewidth=LineSize)
            elif dcvalue == 5:
                ax4.plot(Total_Area.index, Ratio, color=blue_col, linestyle='-.', linewidth=LineSize)
            ax4.set_xlim(t0, tf)
            ax4.set_ylim(0.0, 1.0)
            ax4.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax4.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax4.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax4.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax4.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig4.autofmt_xdate()

    # ------------------------------------------------------------------
    # Build the final legend after all loops have finished
    # ------------------------------------------------------------------
    # Retrieve the automatically labelled C curves and the radar curve
    handles, labels = ax.get_legend_handles_labels()

    handles_C = []
    labels_C = []
    radar_handle = None
    obs_handle = None

    for handle, label in zip(handles, labels):
        if label.startswith('C ='):
            handles_C.append(handle)
            labels_C.append(label)
        elif label == 'Radar':
            radar_handle = handle
        elif label == 'Obs':
            obs_handle = handle

    # Proxy artists describing the different dc values
    handles_dc = [
        Line2D(
            [0], [0],
            color='k',
            linestyle='-',
            linewidth=LineSize,
            label='dc = 2'
        ),
        Line2D(
            [0], [0],
            color='k',
            linestyle='--',
            linewidth=LineSize,
            label='dc = 3'
        ),
        Line2D(
            [0], [0],
            color='k',
            linestyle=':',
            linewidth=LineSize,
            label='dc = 4'
        ),
        Line2D(
            [0], [0],
            color='k',
            linestyle='-.',
            linewidth=LineSize,
            label='dc = 5'
        )
    ]

    labels_dc = ['dc = 2m', 'dc = 3m', 'dc = 4m', 'dc = 5m']

    # Desired legend order:
    # C values, dc values, then radar
    handles_ordered = handles_C + handles_dc
    labels_ordered = labels_C + labels_dc

    if radar_handle is not None:
        handles_ordered.append(radar_handle)
        labels_ordered.append('Radar')
    if obs_handle is not None:
        handles_ordered.append(obs_handle)
        labels_ordered.append('Obs')

    ax.legend(
        handles_ordered,
        labels_ordered,
        fontsize=22
    )

    ax.legend(handles_ordered, labels_ordered, fontsize=22)
    ax2.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower left')
    ax3.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower left')
    ax4.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower left')



    plt.show()