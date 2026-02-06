"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean ub and grounded area over time for chosen C values with both FreeSurface solved
but without flux on lateral BC (final setup for production) to compare results for geometry of 2025, 2011 and 1993

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
    Cvalues= [0.3, 0.4, 0.5]
    # Get 5 shades from the "Blues" colormap (light to dark)
    colors = cmx.get_cmap('viridis_r', len(Cvalues))
    blues = [colors(i) for i in range(len(Cvalues))]
    ###Cases in terms of geometry for which we want a figure:
    Cases = ['1993', '2011', '2025']
    LineStyles = [':','--','-']
    LineSize = 2.4 ###Default value of linewidth
    ###When FS is solved what value of alpha for stabilization for the upper FS
    Alphas = ['10m1']#['10m1','10m2']
    LineSizes = [2.4]#[2.4, 1.1]
    ### Do we plant to plot extra cases for verification (i.e. cases with restart, comparison LUA vs USF_Overload) ?
    VERIF = True
    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=20)
    ax.set_xlabel(r'Time [days]', fontsize=20)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ax.set_title('Mean basal velocity [m/day]', fontsize=21, weight='bold')
    ###~~~~~~~FIG 2 : grounded fraction  ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig2, ax2 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax2.set_ylabel(r'Grounded Area', fontsize=20)
    ax2.set_yscale("linear")
    ax2.set_xlabel(r'Time [days]', fontsize=20)
    ax2.tick_params(labelsize=16)  # fontsize of the tick labels
    ax2.grid(True)
    ax2.grid(alpha=0.5)
    ax2.set_title('Fraction of grounded area', fontsize=21, weight='bold')

    ### START OF SCRIPT###
    ###Loop on each C considered
    for i, (C, blue_col) in enumerate(zip(Cvalues, blues)):
        print('Dealing with C value :', C)
        str_Cvalue = f"{int(C * 10):02d}"
        Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
        ###Loop on year considered
        for j, (Case, LiStyle) in enumerate(zip(Cases, LineStyles)):
            print('Dealing with year :', Case)
            for k, (Alpha, LineSize) in enumerate(zip(Alphas, LineSizes)):
                print('Dealing with Alpha :', Alpha)
                if Case=='2025':
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/3_PRONO_FULLLOAD_FS_NoGMLat/GlacierOut/.')
                    Filename = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp001d_StabFabFSTop_Alpha{}_FSBotStabStdy__bed.dat'.format(str_Cvalue,Alpha)
                    # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
                    pattern = re.compile(r"MyBirch_RS_C{}_BothFS_NoGMLat_tsp001d_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(str_Cvalue))
                    ###Data for control runs
                    Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/2_CONTROL_FS_NoGMLat/GlacierOut/.')
                    Filename_CTRL = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp01d_StabFabFSTop_Alpha{}_FSBotStabStdy__bed.dat'.format(str_Cvalue,Alpha)

                    ### FOR VERIF USF_Overload
                    ###TO BE REMOVED AFTER VERIF
                    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if (C == 0.4 and VERIF):
                        Pathroot_VERIF = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
                        Filename_VERIF = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp001d_CTRL__bed.dat'.format(str_Cvalue)
                        ###load file corresponding to first simu as dataframe
                        print('Opening file:', Filename)
                        Df_wholebed_VERIF = pd.read_csv(Pathroot_VERIF.joinpath(Filename_VERIF), names=Col_Name, sep=r'\s+')
                        ###END OF TO BE REMOVED AFTER VERIF
                        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                else:
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/3_PRONO_FULLLOAD_FS_NoGMLat_{}/GlacierOut/.'.format(Case))
                    if C ==0.5:
                        Filename = 'MyBirch{}_RS_C{}_BothFS_NoGMLat_tsp001d_NoWallTime__bed.dat'.format(Case,str_Cvalue)
                        ###TO BE REMOVED AFTER VERIF
                        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
                        if VERIF:
                            Filename_WithRstt = 'MyBirch{}_RS_C{}_BothFS_NoGMLat_tsp001d_StabFabFSTop_Alpha{}_FSBotStabStdy__bed.dat'.format(Case,str_Cvalue,Alpha)
                            pattern = re.compile(r"MyBirch{}_RS_C{}_BothFS_NoGMLat_tsp001d_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(Case,str_Cvalue))
                        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    else:
                        Filename = 'MyBirch{}_RS_C{}_BothFS_NoGMLat_tsp001d_StabFabFSTop_Alpha{}_FSBotStabStdy__bed.dat'.format(Case,str_Cvalue,Alpha)
                    ###Data for control runs
                    Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/2_CONTROL_FS_NoGMLat_{}/GlacierOut/.'.format(Case))
                    Filename_CTRL = 'MyBirch{}_RS_C{}_BothFS_NoGMLat_tsp01d_StabFabFSTop_Alpha{}_FSBotStabStdy__bed.dat'.format(Case,str_Cvalue,Alpha)
                ###Time step size for prono simu with overload:
                Timestep_size = 0.01 / 365  ### In years !!
                ###load file corresponding to first simu as dataframe
                print('Opening file:', Filename)
                Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')

                ###TO BE REMOVED AFTER VERIF
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if ((Case=='2011' or Case=='1993') and C==0.5 and VERIF):
                    print('Opening file:', Filename_WithRstt)
                    Df_wholebed_WithRstt = pd.read_csv(Pathroot.joinpath(Filename_WithRstt), names=Col_Name, sep=r'\s+')
                    ###Deal with restarts!
                    # Collect restart DataFrames into a list
                    Restart_dfs = []
                    # Gather all restart files and sort them by Rstt (restart number)
                    matches = []
                    for fname in os.listdir(Pathroot):
                        match = pattern.match(fname)
                        if match:
                            Rstt, RsttPos = map(int, match.groups())
                            matches.append((Rstt, RsttPos, fname))
                    matches.sort(key=lambda x: x[0])  # ensure correct order : Restart 1, then restart 2 and so on
                    # Loop efficiently over sorted restart files
                    last_timestep = Df_wholebed_WithRstt["Timestep"].max()
                    last_count = Df_wholebed_WithRstt["Count"].max()
                    ###Now loop on each restart output for current C and year
                    for Rstt, RsttPos, fname in matches:
                        print('Opening file:', fname)
                        df_rstt = pd.read_csv(Pathroot.joinpath(fname), names=Col_Name, sep=r'\s+')
                        ###NOW TREAT THE DATAFRAME OF THE RESTART
                        # Skip first timestep and shift
                        df_rstt = df_rstt[df_rstt["Timestep"] != 1].copy()
                        df_rstt["Timestep"] += last_timestep - 1
                        df_rstt["Count"] += last_count - 1
                        # Update counters
                        last_timestep = df_rstt["Timestep"].max()
                        last_count = df_rstt["Count"].max()
                        ## Add to the full dataframe (first as a list)
                        Restart_dfs.append(df_rstt)

                    # Concatenate once at the end (much faster)
                    if Restart_dfs:  ##Only if Restart_dfs is not empty (there are restarts !)
                        Df_wholebed_WithRstt = pd.concat([Df_wholebed_WithRstt] + Restart_dfs, ignore_index=True)
                ###END OF TO BE REMOVED AFTER VERIF
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                ####NOW LOAD DATA CORRESPONDING TO CONTROL RUNS -> NO RESTART REQUIRED
                ###Time step size for control run:
                Timestep_size_CTRL = 0.1 / 365  ### In years !!
                ###load file as dataframe
                print('Opening file CTRL:', Filename_CTRL)
                Df_wholebed_CTRL = pd.read_csv(Pathroot_CTRL.joinpath(Filename_CTRL), names=Col_Name, sep=r'\s+')
                ###########################################################################
                ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED CASE####
                ###Create a column with time (in years !)
                Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
                Df_wholebed_CTRL['Time'] = Timestep_size_CTRL * Df_wholebed_CTRL['Timestep']
                ###Calculate ub from velocity and normal vector
                ###Start by calculating un (should be zero but there is some residual)
                Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
                Df_wholebed_CTRL['un'] = Df_wholebed_CTRL['u'] * Df_wholebed_CTRL['nx'] + Df_wholebed_CTRL['v'] * Df_wholebed_CTRL['ny'] + Df_wholebed_CTRL['w'] * Df_wholebed_CTRL['nz']
                ###Deduce ub from u and un
                Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
                Df_wholebed_CTRL['ub'] = np.sqrt((Df_wholebed_CTRL['u'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['nx']) ** 2 + (Df_wholebed_CTRL['v'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['ny']) ** 2 + (Df_wholebed_CTRL['w'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['nz']) ** 2)
                ### Calculate mean ub over whole glacier for each time step
                mean_ub_by_time = Df_wholebed.groupby('Time')['ub'].mean()
                mean_ub_by_time_CTRL = Df_wholebed_CTRL.groupby('Time')['ub'].mean()
                ### Show mean ub in Fig 1
                if LiStyle == '-' and k==0:
                    ax.plot(mean_ub_by_time.index*365, mean_ub_by_time/365, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=LineSize)
                else:
                    ax.plot(mean_ub_by_time.index * 365, mean_ub_by_time / 365, color=blue_col, linestyle=LiStyle, linewidth=LineSize)
                ax.plot(mean_ub_by_time_CTRL.index * 365, mean_ub_by_time_CTRL / 365, color=blue_col, linestyle=LiStyle, linewidth=LineSize*0.33)
                ###Dummy plots for legend regarding lateral BC
                if i==0:
                    ax.plot(np.nan,np.nan,label='{}'.format(Case),color='k',linestyle=LiStyle, linewidth=LineSize)
                    if j==0 and k==0:
                        ax.plot(np.nan,np.nan,label='CTRL',color='k',linestyle='-', linewidth=LineSize*0.33)
                ax.set_ylim(0.0, 40.0)
                ax.set_xlim(0, 24.01)
                # Major ticks at specific values
                ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
                # Minor ticks every 1 unit
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)
                ###Plot evolution of grounded area over time in Fig2
                ### Calculate total area over the whole glacier for each time step
                Total_Area = Df_wholebed.groupby('Time')['nodearea'].sum()
                Total_Area_CTRL = Df_wholebed_CTRL.groupby('Time')['nodearea'].sum()
                ### Grounded area
                Total_GroundedArea = (Df_wholebed[Df_wholebed['GM'] == 1].groupby('Time')['nodearea'].sum())
                Total_GroundedArea_CTRL = (Df_wholebed_CTRL[Df_wholebed_CTRL['GM'] == 1].groupby('Time')['nodearea'].sum())
                ### Fraction of grounded area over total area
                Ratio = Total_GroundedArea/Total_Area
                Ratio_CTRL = Total_GroundedArea_CTRL/Total_Area_CTRL
                ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
                # Copy ratio but mask values where Ratio == 1
                Ratio_masked = Ratio.copy()
                Ratio_masked[Ratio_masked == 1] = np.nan
                Ratio_masked_CTRL = Ratio_CTRL.copy()
                Ratio_masked_CTRL[Ratio_masked_CTRL == 1] = np.nan
                ### Plot ratio
                if LiStyle == '-':
                    ax2.plot(Total_Area.index*365, Ratio_masked, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=LineSize)
                else:
                    ax2.plot(Total_Area.index*365, Ratio_masked, color=blue_col,linestyle=LiStyle, linewidth=LineSize)
                ax2.plot(Total_Area_CTRL.index*365, Ratio_masked_CTRL, color=blue_col,linestyle=LiStyle, linewidth=LineSize*0.33)
                ###Dummy plots for legend regarding Case
                if i==0:
                    ax2.plot(np.nan,np.nan,label='Year = {}'.format(Case),color='k',linestyle=LiStyle, linewidth=LineSize)
                    if j==0 and k==0:
                        ax2.plot(np.nan,np.nan,label='CTRL',color='k',linestyle='-', linewidth=LineSize*0.33)
                ax2.set_xlim(0, 24.01)
                # Major ticks at specific values
                ax2.set_xticks([0, 4, 8, 12, 16, 20, 24])
                # Minor ticks every 1 unit
                ax2.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax2.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ####TO REMOVE AFTER VERIF###
                #####To check wether applying USF-Overload instead of LUA function in the Upper BC changes thing
                if VERIF:
                    if Case == '2025' and C == 0.4:
                        ###Create a column with time (in years !)
                        Df_wholebed_VERIF['Time'] = Timestep_size * Df_wholebed_VERIF['Timestep']
                        ###Calculate ub from velocity and normal vector
                        ###Start by calculating un (should be zero but there is some residual)
                        Df_wholebed_VERIF['un'] = Df_wholebed_VERIF['u'] * Df_wholebed_VERIF['nx'] + Df_wholebed_VERIF['v'] * Df_wholebed_VERIF['ny'] + Df_wholebed_VERIF['w'] * Df_wholebed_VERIF['nz']
                        ###Deduce ub from u and un
                        Df_wholebed_VERIF['ub'] = np.sqrt((Df_wholebed_VERIF['u'] - Df_wholebed_VERIF['un'] * Df_wholebed_VERIF['nx']) ** 2 + (Df_wholebed_VERIF['v'] - Df_wholebed_VERIF['un'] * Df_wholebed_VERIF['ny']) ** 2 + (Df_wholebed_VERIF['w'] - Df_wholebed_VERIF['un'] * Df_wholebed_VERIF['nz']) ** 2)
                        ### Calculate mean ub over whole glacier for each time step
                        mean_ub_by_time_VERIF = Df_wholebed_VERIF.groupby('Time')['ub'].mean()
                        ax.plot(mean_ub_by_time_VERIF.index * 365, mean_ub_by_time_VERIF / 365, color='k',linestyle='--', linewidth=LineSize * 0.75)
                        ###Plot evolution of grounded area over time in Fig2
                        ### Calculate total area over the whole glacier for each time step
                        Total_Area_VERIF = Df_wholebed_VERIF.groupby('Time')['nodearea'].sum()
                        ### Grounded area
                        Total_GroundedArea_VERIF = (Df_wholebed_VERIF[Df_wholebed_VERIF['GM'] == 1].groupby('Time')['nodearea'].sum())
                        ### Fraction of grounded area over total area
                        Ratio_VERIF = Total_GroundedArea_VERIF / Total_Area_VERIF
                        Ratio_masked_VERIF = Ratio_VERIF.copy()
                        Ratio_masked_VERIF[Ratio_masked_VERIF == 1] = np.nan
                        ax2.plot(Total_Area_VERIF.index * 365, Ratio_masked_VERIF, color='k', linestyle='--',linewidth=LineSize *0.75)

                    #### To check the effect of having restarts on the solutions
                    if ((Case=='2011' or Case=='1993') and C==0.5):
                        ###Create a column with time (in years !)
                        Df_wholebed_WithRstt['Time'] = Timestep_size * Df_wholebed_WithRstt['Timestep']
                        ###Calculate ub from velocity and normal vector
                        ###Start by calculating un (should be zero but there is some residual)
                        Df_wholebed_WithRstt['un'] = Df_wholebed_WithRstt['u'] * Df_wholebed_WithRstt['nx'] + Df_wholebed_WithRstt['v'] * Df_wholebed_WithRstt['ny'] + Df_wholebed_WithRstt['w'] * Df_wholebed_WithRstt['nz']
                        ###Deduce ub from u and un
                        Df_wholebed_WithRstt['ub'] = np.sqrt((Df_wholebed_WithRstt['u'] - Df_wholebed_WithRstt['un'] * Df_wholebed_WithRstt['nx']) ** 2 + (Df_wholebed_WithRstt['v'] - Df_wholebed_WithRstt['un'] * Df_wholebed_WithRstt['ny']) ** 2 + (Df_wholebed_WithRstt['w'] - Df_wholebed_WithRstt['un'] * Df_wholebed_WithRstt['nz']) ** 2)
                        ### Calculate mean ub over whole glacier for each time step
                        mean_ub_by_time_WithRstt = Df_wholebed_WithRstt.groupby('Time')['ub'].mean()
                        ax.plot(mean_ub_by_time_WithRstt.index * 365, mean_ub_by_time_WithRstt / 365, color='r',linestyle='--', linewidth=LineSize * 0.75)
                        ###Plot evolution of grounded area over time in Fig2
                        ### Calculate total area over the whole glacier for each time step
                        Total_Area_WithRstt = Df_wholebed_WithRstt.groupby('Time')['nodearea'].sum()
                        ### Grounded area
                        Total_GroundedArea_WithRstt = (Df_wholebed_WithRstt[Df_wholebed_WithRstt['GM'] == 1].groupby('Time')['nodearea'].sum())
                        ### Fraction of grounded area over total area
                        Ratio_WithRstt = Total_GroundedArea_WithRstt / Total_Area_WithRstt
                        Ratio_masked_WithRstt = Ratio_WithRstt.copy()
                        Ratio_masked_WithRstt[Ratio_masked_WithRstt == 1] = np.nan
                        ax2.plot(Total_Area_WithRstt.index * 365, Ratio_masked_WithRstt, color='r', linestyle='--',linewidth=LineSize *0.75)



    # --- BUILD FINAL LEGEND WITH CORRECT ORDER ---
    handles, labels = ax.get_legend_handles_labels()
    # Separate C legends (colored solid lines) from BC legends (black dummy)
    handles_C = []
    labels_C = []
    handles_BC = []
    labels_BC = []
    handles_other = []
    labels_other = []

    for h, l in zip(handles, labels):
        if l.startswith("C ="):
            handles_C.append(h)
            labels_C.append(l)
        elif l.startswith("1") or l.startswith("2"):
            handles_BC.append(h)
            labels_BC.append(l)
        else:
            handles_other.append(h)
            labels_other.append(l)

    # Combine them: C first, BC after
    handles_ordered = handles_C + handles_BC + handles_other
    labels_ordered = labels_C + labels_BC + labels_other

    ax.legend(handles_ordered, labels_ordered, fontsize=22)
    ax2.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower right')



    plt.show()