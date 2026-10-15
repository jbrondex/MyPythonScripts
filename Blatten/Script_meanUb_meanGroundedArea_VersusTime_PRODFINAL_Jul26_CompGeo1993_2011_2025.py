"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean velocities (at surface, bed and front) and grounded area over time for a chosen C value and a chosen geometry with both FreeSurface solved
but without flux on lateral BC (final setup for production) and with progressive loading
loading stops either on 19 or 20 or 21 or 22 or 23 of May 25

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
import matplotlib.dates as mdates
from datetime import datetime ###To handle dates

from matplotlib import ticker
from scipy.interpolate import griddata
from matplotlib.ticker import MultipleLocator

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
    Cvalue= 0.3
    dcvalue = 5.0
    Years = ['1993', '2011', '2025']
    ColorBlindFriendly_List = ['k', 'k', 'k'] ##'#984ea3','#999999', '#e41a1c', '#dede00']
    LineSizes = [2.4, 2.4, 2.4]
    LineStyles = [':','--','-']

    ###TEMPORARY TO REMOVE LATER
    IsOld = True ###Use old simu result that stopped randomly for 1993

    StartSimuDay = 16 ####Day of May 2025 corresponding to simulation start
    SimuDuration = 15 ###Simu duration in days
    # StartSimuDay = 19 ####Day of May 2025 corresponding to simulation start
    # SimuDuration = 12 ###Simu duration in days
    EndOfLoadingDays = [23]
    ###Which point of timelapse velocity series do we want to plot in Fig 5 ?
    Points_Timelapse_Ids = [3] ###Ids of points of TimeLapse reconstructed velocities that we want to plot (list with values among 1, 2, 3 and/or 4)
    valid_ids = {1, 2, 3, 4}
    Colors_Points_Timelapse = ['#cc00cc']###List of colors for the points of timelapse worh which we want to have a plot
    ####Rearrangement of parameters to get corresponding file names
    s = f"{Cvalue:.2f}"
    if len(s.rstrip('0').split('.')[1]) == 1:
        str_Cvalue = f"{int(Cvalue * 10):02d}"
    elif len(s.rstrip('0').split('.')[1]) == 2:
        str_Cvalue = f"{int(Cvalue * 100):03d}"

    if dcvalue == 2:
        str_dcvalue = f""
    elif dcvalue <1:
        str_dcvalue = f"dc{int(dcvalue * 10):02d}m_"
    elif dcvalue.is_integer():
        str_dcvalue = f"dc{int(dcvalue):01d}m_"
    else:
        x, y =str(dcvalue).split(".")
        str_dcvalue = f"dc{x}dot{y}m_"

    t0 = pd.Timestamp("2025-05-{} 00:00".format(StartSimuDay))
    tf = t0 + pd.Timedelta(days=SimuDuration)
    CollapseTime = pd.Timestamp("2025-05-28 15:30")
    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'Normal Stress', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']

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
    Points_Timelapse = [[2630784,1139171], [2630758,1139194], [2630712,1139233], [2630681,1139269]]
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

    ###~~~~~~~FIG 5 : Both Front velocity in interferometric radar direction + Horizontal velocity at selected surface points ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig5, ax5 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax5.set_ylabel(r'Ice velocity [m d$^{-1}$]', fontsize=30)
    ax5.set_yscale("linear")
    ax5.tick_params(labelsize=24)  # fontsize of the tick labels
    ax5.grid(True)
    ax5.grid(alpha=0.5)
    ax5.set_title('C = {} and dc ={}m'.format(Cvalue,dcvalue), fontsize=30, weight='bold')

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ### START OF SCRIPT###
    ###Loop on each C considered
    for i, (Year, LineStyle, LineSize, Color) in enumerate(zip(Years, LineStyles, LineSizes, ColorBlindFriendly_List)):
        ##First Data for control runs
        if Year == '2025':
            Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/2_CONTROL_FS_NoGMLat/GlacierOut/.')
            Filename_CTRL = 'MyBirch_RS_C{}_{}tsp01d_CTRL__bed.dat'.format(str_Cvalue,str_dcvalue)
        else:
            Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/2_CONTROL_FS_NoGMLat_{}/GlacierOut/.'.format(Year))
            Filename_CTRL = 'MyBirch{}_RS_C{}_{}tsp01d_CTRL__bed.dat'.format(Year,str_Cvalue,str_dcvalue)
        ###Time step size for control run:
        Timestep_size_CTRL = 0.1 / 365  ### In years !!
        ###load file as dataframe
        print('Opening file CTRL:', Filename_CTRL)
        Df_wholebed_CTRL = pd.read_csv(Pathroot_CTRL.joinpath(Filename_CTRL), names=Col_Name, sep=r'\s+')
        ###Create a column with time (in years !)
        Df_wholebed_CTRL['Time'] = Timestep_size_CTRL * Df_wholebed_CTRL['Timestep']
        ###Convert absolute time in terms of date
        Df_wholebed_CTRL['Date'] = t0 + pd.to_timedelta(Df_wholebed_CTRL['Time'] * 365, unit='D')
        ###Create a dedicated dataframe for each boundary of interest
        Df_wholebed_CTRL_front = Df_wholebed_CTRL[Df_wholebed_CTRL['BC'] == 2]
        Df_wholebed_CTRL_bed = Df_wholebed_CTRL[Df_wholebed_CTRL['BC'] == 3]
        Df_wholebed_CTRL_top = Df_wholebed_CTRL[Df_wholebed_CTRL['BC'] == 4]

        ####MEAN BASAL VELOCITY ####
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed_CTRL_bed['un'] = Df_wholebed_CTRL_bed['u'] * Df_wholebed_CTRL_bed['nx'] + Df_wholebed_CTRL_bed['v'] * Df_wholebed_CTRL_bed['ny'] + Df_wholebed_CTRL_bed['w'] * Df_wholebed_CTRL_bed['nz']
        ###Deduce ub from u and un
        Df_wholebed_CTRL_bed['ub'] = np.sqrt((Df_wholebed_CTRL_bed['u'] - Df_wholebed_CTRL_bed['un'] * Df_wholebed_CTRL_bed['nx']) ** 2 + (Df_wholebed_CTRL_bed['v'] - Df_wholebed_CTRL_bed['un'] * Df_wholebed_CTRL_bed['ny']) ** 2 + (Df_wholebed_CTRL_bed['w'] - Df_wholebed_CTRL_bed['un'] * Df_wholebed_CTRL_bed['nz']) ** 2)
        ###Now calculate the mean basal velocity
        mean_ub_by_time_CTRL = Df_wholebed_CTRL_bed.groupby('Date')['ub'].mean()

        ###GROUNDED AREA FRACTION
        ### Calculate total area over the whole glacier for each time step
        Total_Area_CTRL = Df_wholebed_CTRL_bed.groupby('Date')['nodearea'].sum()
        ### Grounded area
        Total_GroundedArea_CTRL = (Df_wholebed_CTRL_bed[Df_wholebed_CTRL_bed['GM'] == 1].groupby('Date')['nodearea'].sum())
        ### Fraction of grounded area over total area
        Ratio_CTRL = Total_GroundedArea_CTRL / Total_Area_CTRL

        ####FRONT VELOCITY IN RADAR DIRECTION ####
        ###for each node of front, we calculate the vector from radar to node:
        dx = Df_wholebed_CTRL_front['X']-RadarPosition[0]
        dy = Df_wholebed_CTRL_front['Y']-RadarPosition[1]
        dz = Df_wholebed_CTRL_front['Z']-RadarPosition[2]
        ###for each node of front, we calculate norm of this vector:
        r_norm = np.sqrt(dx**2 + dy**2 + dz**2)
        ###Now do the projection through scalar product with unit vector
        Df_wholebed_CTRL_front['uf_radarproj'] = (Df_wholebed_CTRL_front['u']*dx + Df_wholebed_CTRL_front['v']*dy + Df_wholebed_CTRL_front['w']*dz)/r_norm
        ###Now calculate the mean of all nodes at front for each Date
        mean_uf_radarproj_by_time_CTRL = Df_wholebed_CTRL_front.groupby('Date')['uf_radarproj'].mean()

        ####SURFACE VELOCITY####
        ###first calculate mean surface velocity
        Df_wholebed_CTRL_top['us_norm'] = np.sqrt(Df_wholebed_CTRL_top['u'] ** 2 + Df_wholebed_CTRL_top['v'] ** 2 + Df_wholebed_CTRL_top['w'] ** 2)
        ###Now calculate the mean norm of surface velocity for each Date
        mean_us_norm_by_time_CTRL = Df_wholebed_CTRL_top.groupby('Date')['us_norm'].mean()
        ###Now I interpolate the modeled horizontal velocity at the four points of the timelapse camera for each date
        results_CTRL = []
        ##first I will do that for each date available
        for t, group in Df_wholebed_CTRL_top.groupby('Date'):
            ###get x and y of all nodes
            nodes = group[['X', 'Y']].values
            # get horizontal velocities at nodes
            u_vals = group['u'].values
            v_vals = group['v'].values
            # do the interpolation at each of the four points of the timelapse camera
            for id,Point in enumerate(Points_Timelapse):
                u_interp = griddata(nodes, u_vals, Point, method='linear')[0]
                v_interp = griddata(nodes, v_vals, Point, method='linear')[0]
                results_CTRL.append({
                    'Date': t,
                    'Point_id' : id+1, ##id starts at 1
                    'Xpt' : Point[0],
                    'Ypt' : Point[1],
                    'u_interp': u_interp,
                    'v_interp': v_interp
                })
        Df_ModelInterpAtTimeLapsePoints_CTRL=pd.DataFrame(results_CTRL)
        ###Calculate norm of horizontal velocity
        Df_ModelInterpAtTimeLapsePoints_CTRL['u_interp_norm']=np.sqrt(Df_ModelInterpAtTimeLapsePoints_CTRL['u_interp'] ** 2 + Df_ModelInterpAtTimeLapsePoints_CTRL['v_interp'] ** 2)

        ####NOW LOAD DATA CORRESPONDING TO PRONOSTIC RUNS WITH OVERLOAD
        ###For now I hard code the line below, to change later if needed
        EndDay = EndOfLoadingDays[0]
        print('Dealing with simu for year ', Year, 'with overload finished on :', EndDay, 'of May')
        if Year == '2025':
            Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
            Filename = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad{}May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue,EndDay)
        elif (Year == '1993' and IsOld):
            Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat_{}/GlacierOut/.'.format(Year))
            Filename = 'MyBirch{}_RS_C{}_{}tsp001d_Start16May_StopLoad{}May_Stop31May__bed_OLD.dat'.format(Year, str_Cvalue,str_dcvalue,EndDay)
        else:
            Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat_{}/GlacierOut/.'.format(Year))
            Filename = 'MyBirch{}_RS_C{}_{}tsp001d_Start16May_StopLoad{}May_Stop31May__bed.dat'.format(Year, str_Cvalue,str_dcvalue,EndDay)
        ###Time step size for prono simu with overload:
        Timestep_size = 0.01 / 365  ### In years !!
        ###load file corresponding to first simu as dataframe
        print('Opening file:', Filename)
        Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
        ###########################################################################
        ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED END OF OVERLOAD DAY
        ###Create a column with time (in years !)
        Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
        ###Convert absolute time in terms of date
        Df_wholebed['Date'] = t0 + pd.to_timedelta(Df_wholebed['Time'] * 365, unit='D')
        ###Create a dedicated dataframe for each boundary of interest
        Df_wholebed_front = Df_wholebed[Df_wholebed['BC'] == 2]
        Df_wholebed_bed = Df_wholebed[Df_wholebed['BC'] == 3]
        Df_wholebed_top = Df_wholebed[Df_wholebed['BC'] == 4]

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
        ###first calculate mean surface velo
        Df_wholebed_top['us_norm'] = np.sqrt(Df_wholebed_top['u'] ** 2 + Df_wholebed_top['v'] ** 2 + Df_wholebed_top['w'] ** 2)
        ###Now calculate the mean norm of surface velocity for each Date
        mean_us_norm_by_time = Df_wholebed_top.groupby('Date')['us_norm'].mean()
        ###Now I interpolate the modeled horizontal velocity at the four points of the timelapse camera for each date
        results = []
        #first I will do that for each date available
        for t, group in Df_wholebed_top.groupby('Date'):
            ###get x and y of all nodes
            nodes = group[['X', 'Y']].values
            # get horizontal velocities at nodes
            u_vals = group['u'].values
            v_vals = group['v'].values
            # do the interpolation at each of the four points of the timelapse camera
            for id, Point in enumerate(Points_Timelapse):
                u_interp = griddata(nodes, u_vals, Point, method='linear')[0]
                v_interp = griddata(nodes, v_vals, Point, method='linear')[0]
                results.append({
                    'Date': t,
                    'Point_id': id+1, ##id starts at 1
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
        ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, label='{}'.format(Year), color=Color, linestyle=LineStyle, linewidth=LineSize)
        if LineStyle == '-':
            ax.plot(mean_ub_by_time_CTRL.index, mean_ub_by_time_CTRL / 365, label='CTRL', color='0.3', linestyle=LineStyle, linewidth=LineSize*0.4)
            ax.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs() / 1000, color='b', label="Radar obs. May 25 (front)", linestyle='--',linewidth=LineSize * 1.6)
            ###Add vertical line for collapse day and annotate
            ax.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
            ax.text(CollapseTime, 40, 'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
        else:
            ax.plot(mean_ub_by_time_CTRL.index, mean_ub_by_time_CTRL / 365, color='0.3', linestyle=LineStyle, linewidth=LineSize*0.4)
        ###Dummy plots for legend regarding lateral B
        ax.set_ylim(0.0, 40.0)
        ax.set_xlim(t0, tf)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
        # Add grid for both major and minor ticks
        ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
        fig1.autofmt_xdate()



        ###Fig 3: Velocity at front in radar direction
        ax3.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs()/365, label= '{}'.format(Year), color=Color,linestyle=LineStyle, linewidth=LineSize)
        if LineStyle == '-':
            ax3.plot(mean_uf_radarproj_by_time_CTRL.index, mean_uf_radarproj_by_time_CTRL.abs()/365, label= 'CTRL', color='0.3',linestyle=LineStyle, linewidth=LineSize*0.4)
            ax3.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs() / 1000, color='b', label="Radar obs. May 25 (front)", linestyle='--',linewidth=LineSize * 1.6)
            ###Add vertical line for collapse day and annotate
            ax3.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
            ax3.text(CollapseTime, 40,'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
        else:
            ax3.plot(mean_uf_radarproj_by_time_CTRL.index, mean_uf_radarproj_by_time_CTRL.abs()/365, color='0.3',linestyle=LineStyle, linewidth=LineSize*0.4)
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
        ax4.plot(Total_Area.index, Ratio, label='{}'.format(Year), color=Color,linestyle=LineStyle, linewidth=LineSize)
        if LineStyle == '-':
            ax4.plot(Total_Area_CTRL.index, Ratio_CTRL, label='CTRL', color='0.3', linestyle=LineStyle, linewidth=LineSize*0.4)
            ###Add vertical line for collapse day
            ax4.axvline(CollapseTime,linestyle='--',linewidth=2.7, color= 'r')
            ax4.text(CollapseTime, 40, 'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
        else:
            ax4.plot(Total_Area_CTRL.index, Ratio_CTRL, color='0.3', linestyle=LineStyle, linewidth=LineSize*0.4)
        ax4.set_xlim(t0, tf)
        ax4.set_ylim(0.0, 1.0)
        ax4.xaxis.set_major_locator(mdates.DayLocator(interval=2))
        ax4.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax4.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
        # Add grid for both major and minor ticks
        ax4.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
        ax4.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
        fig4.autofmt_xdate()

        ###Fig 5: Velocity at front in radar direction + Horizontal velocity at selected surface points + mean ub ~~~~~~~~~~~~~~~
        ### TimeLapse Velocities Modeled versus Observed
        for Points_Timelapse_Id, colstar in zip(Points_Timelapse_Ids, Colors_Points_Timelapse):
            if len(Points_Timelapse_Ids) != len(Colors_Points_Timelapse):
                raise ValueError(
                    "Points_Timelapse_Ids and Colors_Points_Timelapse "
                    "must have the same length."
                )
            if Points_Timelapse_Id not in valid_ids:
                raise ValueError(f"Invalid value in Points_Timelapse_Ids: '{Points_Timelapse_Id}'. "
                                 "There are 4 points at which we have data so should be 1, 2, 3 or 4.")
            Point = Points_Timelapse[Points_Timelapse_Id - 1]
            subset2 = Df_ModelInterpAtTimeLapsePoints[Df_ModelInterpAtTimeLapsePoints['Point_id'] == Points_Timelapse_Id]
            ax5.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365, label= '{}'.format(Year), color=Color,linestyle=LineStyle, linewidth=LineSize)
            ax5.plot(subset2['Date'], subset2['u_interp_norm'] / 365, color=colstar, linestyle=LineStyle, linewidth=LineSize)
            if LineStyle == '-':
                ax5.plot(subset2['Date'], subset2['u_interp_norm'] / 365,label="Modelled velocity (surface point)", color=colstar, linestyle=LineStyle,linewidth=LineSize)
                subset_Obs = Df_Timelapse[(Df_Timelapse['X0'] == Point[0]) & (Df_Timelapse['Y0'] == Point[1])]
                ax5.scatter(subset_Obs['t_mid'], subset_Obs['horizontal_velo'], label="Time-lapse obs May 25 (surface point)", color=colstar, s=20)
                ax5.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs() / 365, label="Modelled front velocity", color=Color, linestyle=LineStyle, linewidth=LineSize)
                ax5.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs() / 1000, color='b',label="Radar obs. May 25 (front)", linestyle='--', linewidth=LineSize * 1.6)
                ###Add vertical line for collapse day and annotate
                ax5.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
                ax5.text(CollapseTime, 40,'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')

            ###Dummy plots for legend regarding lateral B
            ax5.set_ylim(0.0, 100.0)
            ax5.set_xlim(t0, tf)
            ax5.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax5.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax5.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax5.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax5.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig5.autofmt_xdate()

    # --- BUILD FINAL LEGEND WITH CORRECT ORDER ---
    handles, labels = ax.get_legend_handles_labels()
    # Separate C legends (colored solid lines) from BC legends (black dummy)
    handles_Yr = []
    labels_Yr = []
    handles_CTRL = []
    labels_CTRL = []
    handles_Radar = []
    labels_Radar = []
    handles_other = []
    labels_other = []

    for h, l in zip(handles, labels):
        if l.startswith("1") or l.startswith("2"): ###Year
            handles_Yr.append(h)
            labels_Yr.append(l)
        elif l.startswith("CTR"):
            handles_CTRL.append(h)
            labels_CTRL.append(l)
        elif l.startswith("Ra"):
            handles_Radar.append(h)
            labels_Radar.append(l)
        else:
            handles_other.append(h)
            labels_other.append(l)

    # Combine them: C first, BC after
    handles_ordered = handles_Yr + handles_CTRL + handles_Radar + handles_other
    labels_ordered = labels_Yr + labels_CTRL + labels_Radar + labels_other

    legend1 = ax.legend(handles_ordered, labels_ordered,  fontsize=22)
    legend3 = ax3.legend(handles_ordered, labels_ordered,  fontsize=22, loc='upper left')

    handles_ordered4 = handles_Yr + handles_CTRL
    labels_ordered4 = labels_Yr + labels_CTRL
    legend4 = ax4.legend(handles_ordered4, labels_ordered4,  fontsize=22, loc='lower left')

    ####ONE SPECIFIC HANDLING OF LEGEND FOR FIGURE 5
    # Retrieve handles and labels specifically from ax5
    handles5, labels5 = ax5.get_legend_handles_labels()

    main_labels = [
        "Modelled velocity (surface point)",
        "Time-lapse obs May 25 (surface point)",
        "Modelled front velocity",
        "Radar obs. May 25 (front)",

    ]

    handles_main = []
    labels_main = []

    handles_Yr = []
    labels_Yr = []

    # Avoid repeated entries caused by the loop over timelapse points
    seen_main = set()
    seen_Yr = set()

    for handle, label in zip(handles5, labels5):

        # Main observation/model entries
        if label in main_labels and label not in seen_main:
            handles_main.append(handle)
            labels_main.append(label)
            seen_main.add(label)

        # Labels corresponding to the year
        elif (label.startswith("1") or label.startswith("2")) and label not in seen_Yr:
            handles_Yr.append(handle)
            labels_Yr.append(label)
            seen_Yr.add(label)

    # Ligne invisible servant de sous-titre
    dummy = Line2D([], [], linestyle="none")

    handles5_merged = handles_main + [dummy] + handles_Yr
    labels5_merged = labels_main + ["Geometry of Year:"] + labels_Yr
    # First legend: observations and full-load simulation
    legend5_main = ax5.legend(
        handles5_merged,
        labels5_merged,
        fontsize=20,
        loc="upper left"
    )
    # Rendre le handle du sous-titre invisible
    legend5_main .legend_handles[len(handles_main)].set_visible(False)
    # Mettre le sous-titre en gras (optionnel)
    legend5_main.get_texts()[len(handles_main)].set_weight("bold")




    plt.show()