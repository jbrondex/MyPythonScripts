"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean ub and grounded area over time for a chosen C value and a chosen geometry with both FreeSurface solved
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
    str_Cvalue = f"{int(Cvalue * 10):02d}"
    Year = '2025'
    StartSimuDay = 16 ####Day of May 2025 corresponding to simulation start
    SimuDuration = 15 ###Simu duration in days
    # StartSimuDay = 19 ####Day of May 2025 corresponding to simulation start
    # SimuDuration = 12 ###Simu duration in days
    ###Do we want to keep only nodes that belongs to both bed and front boundaries ??? (For comparison to radar observations while using output made at bed only)
    KeepOnlyFrontNodes = True
    EndOfLoadingDays = [19,20,21,22,23]
    t0 = pd.Timestamp("2025-05-{} 00:00".format(StartSimuDay))
    tf = t0 + pd.Timedelta(days=SimuDuration)
    CollapseTime = pd.Timestamp("2025-05-28 15:30")
    #EndOfLoadingDays = [23]
    # colors = ['c'] #['c','r','m']##,'b', 'k']
    ###In this order: blue/orange/green/brown/purple/grey/red/yellow
    ColorBlindFriendly_List = ['#377eb8', '#ff7f00', '#4daf4a', '#a65628', 'k'] ##'#984ea3','#999999', '#e41a1c', '#dede00']
    LineSizes = [1.8, 1.8, 1.8, 1.8, 2.4]
    #ColorBlindFriendly_List = ['k'] ##'#984ea3','#999999', '#e41a1c', '#dede00']
    #LineSizes = [3.4]
    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']


    ###TREAT OBSERVATIONS###
    LastReliableRadarData = pd.Timestamp("2025-05-25 06:00") ###Radar data is not reliable anymore after 6 am on the 25th of May 25
    Col_Name_Radar = ['Time','Fullfront_6h','Middlefront_6h','Upperfront_6h','Fullfront_24h','Middlefront_24h','Upperfront_24h','Fullfront_RWB','Middlefront_RWB','Upperfront_RWB','FullMean']
    Pathroot_Radar = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/Obs/.')
    Filename_Radar = 'glacier_front_RadarData_mmperday.csv'
    ###BE WARE : Radar data are in mm/day
    Df_Radar = pd.read_csv(Pathroot_Radar.joinpath(Filename_Radar), header=0, names=Col_Name_Radar)
    ###Convert time of csv file in panda date
    Df_Radar['Date'] = pd.to_datetime(Df_Radar["Time"], format="%d.%m.%y %H:%M")
    ###Cut all data that are not reliable (i.e. aquired later than last reliable data)
    Df_Radar =  Df_Radar[Df_Radar['Date']<LastReliableRadarData]
    ###Calculate a mean between the three processing method (6h vs 24h vs RWB) for the full front
    Df_Radar["mean_Fullfront"] = Df_Radar[["Fullfront_6h", "Fullfront_24h", "Fullfront_RWB"]].mean(axis=1, skipna=True)

    ###If we want to keep only the nodes belonging to both bed and front boundary, them we need to identify thoses nodes
    if KeepOnlyFrontNodes :
        Col_Name_ForFront = ['Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'zsRock', 'zsIce']
        Pathroot_ForFront = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/1_Mesh/MESH2025/GlacierOut/.')
        Filename_ForFront = 'MyBirch_Init_DEMCorrected_GL_FS_SLFront__bed.dat'
        Df_ForFront = pd.read_csv(Pathroot_ForFront.joinpath(Filename_ForFront), names=Col_Name_ForFront, sep=r'\s+')
        ###Identify nodes of bed which have the same X and Y as nodes of front
        nodes_with_both = (Df_ForFront[Df_ForFront["BC"].isin([2, 3])].groupby(["X", "Y"])["BC"].nunique().eq(2))
        ###Keep only nodes of bed but fo which X,Y correspond to front
        df_bc3 = Df_ForFront[Df_ForFront["BC"] == 3]
        filtered_nodes = df_bc3.set_index(["X", "Y"]).loc[nodes_with_both[nodes_with_both].index].reset_index()
    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    if KeepOnlyFrontNodes :
        ax.set_ylabel(r'$\langle u_\mathrm{f} \rangle$ [m/day]', fontsize=30)
    else:
        ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=30)
    #ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/jour]', fontsize=30)
    ax.tick_params(labelsize=24)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    #ax.set_title('Vitesse basale moyenne [m/jour]', fontsize=30, weight='bold')
    # if KeepOnlyFrontNodes :
    #     ax.set_title('Mean velocity at front [m/day]', fontsize=30, weight='bold')
    # else:
    #     ax.set_title('Mean basal velocity [m/day]', fontsize=30, weight='bold')
    ###~~~~~~~FIG 2 : grounded fraction  ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig2, ax2 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax2.set_ylabel(r'Grounded Area', fontsize=30)
    ax2.set_yscale("linear")
    ax2.tick_params(labelsize=24)  # fontsize of the tick labels
    ax2.grid(True)
    ax2.grid(alpha=0.5)
    ax2.set_title('Fraction of grounded area', fontsize=30, weight='bold')
    ### START OF SCRIPT###
    ###Loop on each C considered
    for i, (EndDay, blue_col, LineSize) in enumerate(zip(EndOfLoadingDays, ColorBlindFriendly_List, LineSizes)):
        ###First Data for control runs
        if i ==0:
            Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/2_CONTROL_FS_NoGMLat/GlacierOut/.')
            Filename_CTRL = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp01d_StabFabFSTop_Alpha10m1_FSBotStabStdy__bed.dat'.format(str_Cvalue)
            ###Time step size for control run:
            Timestep_size_CTRL = 0.1 / 365  ### In years !!
            ###load file as dataframe
            print('Opening file CTRL:', Filename_CTRL)
            Df_wholebed_CTRL = pd.read_csv(Pathroot_CTRL.joinpath(Filename_CTRL), names=Col_Name, sep=r'\s+')
            ###Create a column with time (in years !)
            Df_wholebed_CTRL['Time'] = Timestep_size_CTRL * Df_wholebed_CTRL['Timestep']
            ###Convert absolute time in terms of date
            Df_wholebed_CTRL['Date'] = t0 + pd.to_timedelta(Df_wholebed_CTRL['Time'] * 365, unit='D')
            ###Calculate ub from velocity and normal vector
            ###Start by calculating un (should be zero but there is some residual)
            Df_wholebed_CTRL['un'] = Df_wholebed_CTRL['u'] * Df_wholebed_CTRL['nx'] + Df_wholebed_CTRL['v'] * Df_wholebed_CTRL['ny'] + Df_wholebed_CTRL['w'] * Df_wholebed_CTRL['nz']
            ###Deduce ub from u and un
            Df_wholebed_CTRL['ub'] = np.sqrt((Df_wholebed_CTRL['u'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['nx']) ** 2 + (Df_wholebed_CTRL['v'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['ny']) ** 2 + (Df_wholebed_CTRL['w'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['nz']) ** 2)
            ### Calculate mean ub over whole glacier for each time step
            if KeepOnlyFrontNodes:
                ###filter nodes that belongs to the front (while being also on the bed)
                Df_wholebed_FrontFiltered_CTRL = Df_wholebed_CTRL[Df_wholebed_CTRL["NodeNumber"].isin(filtered_nodes["NodeNumber"])]
                mean_ub_by_time_CTRL = Df_wholebed_FrontFiltered_CTRL.groupby('Date')['ub'].mean()
            else:
                mean_ub_by_time_CTRL = Df_wholebed_CTRL.groupby('Date')['ub'].mean()
            ### Calculate total area over the whole glacier for each time step
            Total_Area_CTRL = Df_wholebed_CTRL.groupby('Date')['nodearea'].sum()
            ### Grounded area
            Total_GroundedArea_CTRL = (Df_wholebed_CTRL[Df_wholebed_CTRL['GM'] == 1].groupby('Date')['nodearea'].sum())
            ### Fraction of grounded area over total area
            Ratio_CTRL = Total_GroundedArea_CTRL / Total_Area_CTRL

        ####NOW LOAD DATA CORRESPONDING TO PRONOSTIC RUNS WITH OVERLOAD
        print('Dealing with simu with overload finished on :', EndDay, 'of May')
        Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/5_PRONO_PROGRESSIVELOAD_FS_NoGMLat_MaxLoadFromDEMs/GlacierOut/.')
#        Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
        Filename = 'MyBirch_RS_C{}_tsp001d_Start{}May_StopLoad{}May_Stop{}May__bed.dat'.format(str_Cvalue,StartSimuDay,EndDay,StartSimuDay+SimuDuration)
        # Filename = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp001d_StopLoad{}May__bed.dat'.format(str_Cvalue,EndDay)
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
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
        ###Deduce ub from u and un
        Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
        ### Calculate mean ub over whole glacier for each time step
        if KeepOnlyFrontNodes :
            ###filter nodes that belongs to the front (while being also on the bed)
            Df_wholebed_FrontFiltered = Df_wholebed[Df_wholebed["NodeNumber"].isin(filtered_nodes["NodeNumber"])]
            mean_ub_by_time = Df_wholebed_FrontFiltered.groupby('Date')['ub'].mean()
        else:
            mean_ub_by_time = Df_wholebed.groupby('Date')['ub'].mean()
        ### Show mean ub in Fig 1
        if blue_col=='k':
            ax.plot(mean_ub_by_time.index, mean_ub_by_time/365, label= 'C = {}'.format(Cvalue), color=blue_col,linestyle='-', linewidth=LineSize)
        else:
            ax.plot(mean_ub_by_time.index, mean_ub_by_time / 365, label='{} May'.format(EndDay), color=blue_col, linestyle='-', linewidth=LineSize)
        if i == 0:
            # ax.plot(mean_ub_by_time_CTRL.index, mean_ub_by_time_CTRL / 365, color='k', label= 'CTRL', linestyle='--', linewidth=LineSize)
            ax.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs()/1000, color='b', label= 'Radar', linestyle='--', linewidth=LineSize*2)
            # ax.plot(Df_Radar['Date'], Df_Radar['FullMean'].abs()/1000, color='b', linestyle='--', linewidth=LineSize*2)
            ###Add vertical line for collapse day and annotate
            ax.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
            ax.text(CollapseTime, 40,'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
            #ax.text(CollapseTime, 40,'Effondrement', fontsize=26, color='r', rotation=90, va='top', ha='right')
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
        ###FOR THE FIGURE OF GROUNDED AREA FRACTION
        ###Plot evolution of grounded area over time in Fig2
        ### Calculate total area over the whole glacier for each time step
        Total_Area = Df_wholebed.groupby('Date')['nodearea'].sum()
        ### Grounded area
        Total_GroundedArea = (Df_wholebed[Df_wholebed['GM'] == 1].groupby('Date')['nodearea'].sum())
        ### Fraction of grounded area over total area
        Ratio = Total_GroundedArea/Total_Area
        ### Plot ratio
        if blue_col=='k':
            ax2.plot(Total_Area.index, Ratio, label= 'C = {}'.format(Cvalue), color=blue_col,linestyle='-', linewidth=LineSize)
        else:
            ax2.plot(Total_Area.index, Ratio, label='{} May'.format(EndDay), color=blue_col, linestyle='-', linewidth=LineSize)
        if i==0:
            ax2.plot(Total_Area_CTRL.index, Ratio_CTRL, color='k', label= 'CTRL', linestyle='--', linewidth=LineSize)
            ###Add vertical line for collapse day
            ax2.axvline(CollapseTime,linestyle='--',linewidth=2.7, color= 'r')
            ax2.text(CollapseTime, 40, 'Collapse', fontsize=26, color='r', rotation=90, va='top', ha='right')
        ax2.set_xlim(t0, tf)
        ax2.set_ylim(0.0, 1.0)
        ax2.xaxis.set_major_locator(mdates.DayLocator(interval=2))
        ax2.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
        # Add grid for both major and minor ticks
        ax2.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
        ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
        fig2.autofmt_xdate()

    # --- BUILD FINAL LEGEND WITH CORRECT ORDER ---
    handles, labels = ax.get_legend_handles_labels()
    # Separate C legends (colored solid lines) from BC legends (black dummy)
    handles_C = []
    labels_C = []
    handles_BC = []
    labels_BC = []
    handles_Radar = []
    labels_Radar = []
    handles_other = []
    labels_other = []

    for h, l in zip(handles, labels):
        if l.startswith("C ="):
            handles_C.append(h)
            labels_C.append(l)
        elif l.startswith("1") or l.startswith("2"):
            handles_BC.append(h)
            labels_BC.append(l)
        elif l.startswith("Ra"):
            handles_Radar.append(h)
            labels_Radar.append(l)
        else:
            handles_other.append(h)
            labels_other.append(l)

    # Combine them: C first, BC after
    handles_ordered = handles_C + handles_BC + handles_Radar + handles_other
    labels_ordered = labels_C + labels_BC + labels_Radar + labels_other
    handles_ordered2 = handles_C + handles_BC + handles_other
    labels_ordered2 = labels_C + labels_BC + labels_other

    # legend1 = ax.legend(handles_ordered, labels_ordered, title="Last day of loading:", fontsize=22)
    # legend2 = ax2.legend(handles_ordered2, labels_ordered2, title="Last day of loading:", fontsize=22, loc='lower left')
    # legend1.get_title().set_fontsize(22)
    # legend2.get_title().set_fontsize(22)

    plt.show()