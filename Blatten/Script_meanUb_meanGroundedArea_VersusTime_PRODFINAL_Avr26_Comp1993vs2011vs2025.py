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
    # Cases = ['1993', '2011', '2025']
    # LineStyles = [':','--','-']
    Cases = ['2025']
    LineStyles = ['--']
    LineSize = 2.4 ###Default value of linewidth
    ###When FS is solved what value of alpha for stabilization for the upper FS
    Alphas = ['10m1']#['10m1','10m2']
    LineSizes = [2.4]#[2.4, 1.1]
    ### Over which boundary do we want to plot the velocities ?
    BoundaryToPlot = 'Bed' ###Can be 'Front', 'Top' or 'Bed'
    ###BE WARE : For now I consider norm of basal (tangential) velocity at the bed, while I consider norm of full veloticy vector
    ### at the front and surface
    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax.set_xlabel(r'Time [days]', fontsize=20)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    if BoundaryToPlot == 'Bed':
        ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=20)
        ax.set_title('Mean basal velocity [m/day]', fontsize=21, weight='bold')
    elif BoundaryToPlot == 'Front':
        ax.set_ylabel(r'$\langle u_\mathrm{f} \rangle$ [m/day]', fontsize=20)
        ax.set_title('Mean velocity at front [m/day]', fontsize=21, weight='bold')
    elif BoundaryToPlot == 'Top':
        ax.set_ylabel(r'$\langle u_\mathrm{s} \rangle$ [m/day]', fontsize=20)
        ax.set_title('Mean surface velocity [m/day]', fontsize=21, weight='bold')
    else:
        raise ValueError(f"Unknown BoundaryToPlot value: {BoundaryToPlot}. Should be 'Front' or 'Bed' or 'Top'")
    ###~~~~~~~FIG 2 : grounded fraction  ~~~~~~~~~~~~~~~
    if BoundaryToPlot == 'Bed':
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
        ####Rearrangement of parameters to get corresponding file names
        s = f"{C:.2f}"
        if len(s.rstrip('0').split('.')[1]) == 1:
            str_Cvalue = f"{int(C * 10):02d}"
        elif len(s.rstrip('0').split('.')[1]) == 2:
            str_Cvalue = f"{int(C * 100):03d}"

        Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
        Col_Name2 = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'Normal Stress', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
        ###Loop on year considered
        for j, (Case, LiStyle) in enumerate(zip(Cases, LineStyles)):
            print('Dealing with year :', Case)
            for k, (Alpha, LineSize) in enumerate(zip(Alphas, LineSizes)):
                if Case=='2025':
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/3_PRONO_FULLLOAD_FS_NoGMLat/GlacierOut/.')
                    Filename = 'MyBirch_RS_C{}_tsp001d_FullLoad__bed.dat'.format(str_Cvalue)
                    ###Data for control runs
                    Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/2_CONTROL_FS_NoGMLat/GlacierOut/.')
                    Filename_CTRL = 'MyBirch_RS_C{}_tsp01d_CTRL__bed.dat'.format(str_Cvalue)
                    Pathroot2 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/3_PRONO_FULLLOAD_FS_NoGMLat/GlacierOut/.')
                    Filename2 = 'MyBirch_RS_C{}_dc2m_tsp001d_FullLoad__bed.dat'.format(str_Cvalue)
                    ###Data for control runs
                    Pathroot_CTRL2 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/2_CONTROL_FS_NoGMLat/GlacierOut/.')
                    Filename_CTRL2 = 'MyBirch_RS_C{}_dc2m_tsp01d_CTRL__bed.dat'.format(str_Cvalue)
                else:
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/3_PRONO_FULLLOAD_FS_NoGMLat_{}/GlacierOut/.'.format(Case))
                    Filename = 'MyBirch{}_RS_C{}_tsp001d_FullLoad__bed.dat'.format(Case,str_Cvalue)
                    ###Data for control runs
                    Pathroot_CTRL = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/2_CONTROL_FS_NoGMLat_{}/GlacierOut/.'.format(Case))
                    Filename_CTRL = 'MyBirch{}_RS_C{}_tsp01d_CTRL__bed.dat'.format(Case,str_Cvalue)
                ###Time step size for prono simu with overload:
                Timestep_size = 0.01 / 365  ### In years !!
                ###load file corresponding to first simu as dataframe
                print('Opening file:', Filename)
                Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
                Df_wholebed2 = pd.read_csv(Pathroot2.joinpath(Filename2), names=Col_Name2, sep=r'\s+')
                ###Keep only boundary of interest
                if BoundaryToPlot == 'Front':
                    Df_wholebed = Df_wholebed[Df_wholebed['BC'] == 2]
                elif BoundaryToPlot == 'Bed':
                    Df_wholebed = Df_wholebed[Df_wholebed['BC'] == 3]
                    Df_wholebed2 = Df_wholebed2[Df_wholebed2['BC'] == 3]
                elif BoundaryToPlot == 'Top' :
                    Df_wholebed = Df_wholebed[Df_wholebed['BC'] == 4]
                else:
                    raise ValueError(f"Unknown BoundaryToPlot value: {BoundaryToPlot}. Should be 'Front' or 'Bed' or 'Top'")

                ####NOW LOAD DATA CORRESPONDING TO CONTROL RUNS -> NO RESTART REQUIRED
                ###Time step size for control run:
                Timestep_size_CTRL = 0.1 / 365  ### In years !!
                ###load file as dataframe
                print('Opening file CTRL:', Filename_CTRL)
                Df_wholebed_CTRL = pd.read_csv(Pathroot_CTRL.joinpath(Filename_CTRL), names=Col_Name, sep=r'\s+')
                Df_wholebed_CTRL2 = pd.read_csv(Pathroot_CTRL2.joinpath(Filename_CTRL2), names=Col_Name2, sep=r'\s+')
                ###Keep only boundary of interest
                if BoundaryToPlot == 'Front':
                    Df_wholebed_CTRL  = Df_wholebed_CTRL [Df_wholebed_CTRL ['BC'] == 2]
                elif BoundaryToPlot == 'Bed':
                    Df_wholebed_CTRL  = Df_wholebed_CTRL [Df_wholebed_CTRL ['BC'] == 3]
                    Df_wholebed_CTRL2  = Df_wholebed_CTRL2 [Df_wholebed_CTRL2 ['BC'] == 3]
                elif BoundaryToPlot == 'Top' :
                    Df_wholebed_CTRL  = Df_wholebed[Df_wholebed_CTRL ['BC'] == 4]
                else:
                    raise ValueError(f"Unknown BoundaryToPlot value: {BoundaryToPlot}. Should be 'Front' or 'Bed' or 'Top'")
                ###########################################################################
                ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED CASE####
                ###Create a column with time (in years !)
                Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
                Df_wholebed_CTRL['Time'] = Timestep_size_CTRL * Df_wholebed_CTRL['Timestep']
                Df_wholebed2['Time'] = Timestep_size * Df_wholebed2['Timestep']
                Df_wholebed_CTRL2['Time'] = Timestep_size_CTRL * Df_wholebed_CTRL2['Timestep']
                if BoundaryToPlot =='Bed':
                    ###Calculate ub from velocity and normal vector
                    ###Start by calculating un (should be zero but there is some residual)
                    Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
                    Df_wholebed_CTRL['un'] = Df_wholebed_CTRL['u'] * Df_wholebed_CTRL['nx'] + Df_wholebed_CTRL['v'] * Df_wholebed_CTRL['ny'] + Df_wholebed_CTRL['w'] * Df_wholebed_CTRL['nz']
                    ###Deduce ub from u and un
                    Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
                    Df_wholebed_CTRL['ub'] = np.sqrt((Df_wholebed_CTRL['u'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['nx']) ** 2 + (Df_wholebed_CTRL['v'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['ny']) ** 2 + (Df_wholebed_CTRL['w'] - Df_wholebed_CTRL['un'] * Df_wholebed_CTRL['nz']) ** 2)
                    ####Samething with new implementation
                    Df_wholebed2['un'] = Df_wholebed2['u'] * Df_wholebed2['nx'] + Df_wholebed2['v'] * Df_wholebed2['ny'] + Df_wholebed2['w'] * Df_wholebed2['nz']
                    Df_wholebed_CTRL2['un'] = Df_wholebed_CTRL2['u'] * Df_wholebed_CTRL2['nx'] + Df_wholebed_CTRL2['v'] * Df_wholebed_CTRL2['ny'] + Df_wholebed_CTRL2['w'] * Df_wholebed_CTRL2['nz']
                    ###Deduce ub from u and un
                    Df_wholebed2['ub'] = np.sqrt((Df_wholebed2['u'] - Df_wholebed2['un'] * Df_wholebed2['nx']) ** 2 + (Df_wholebed2['v'] - Df_wholebed2['un'] * Df_wholebed2['ny']) ** 2 + (Df_wholebed2['w'] - Df_wholebed2['un'] * Df_wholebed2['nz']) ** 2)
                    Df_wholebed_CTRL2['ub'] = np.sqrt((Df_wholebed_CTRL2['u'] - Df_wholebed_CTRL2['un'] * Df_wholebed_CTRL2['nx']) ** 2 + (Df_wholebed_CTRL2['v'] - Df_wholebed_CTRL2['un'] * Df_wholebed_CTRL2['ny']) ** 2 + (Df_wholebed_CTRL2['w'] - Df_wholebed_CTRL2['un'] * Df_wholebed_CTRL2['nz']) ** 2)

                elif BoundaryToPlot =='Front' or BoundaryToPlot =='Top':
                    ###Deduce ub from u and un
                    Df_wholebed['ub'] = np.sqrt(Df_wholebed['u']** 2 + Df_wholebed['v']** 2 + Df_wholebed['w']** 2)
                    Df_wholebed_CTRL['ub'] = np.sqrt(Df_wholebed_CTRL['u'] ** 2 + Df_wholebed_CTRL['v']** 2 + Df_wholebed_CTRL['w']** 2)
                ### Calculate mean ub over whole glacier for each time step
                mean_ub_by_time = Df_wholebed.groupby('Time')['ub'].mean()
                mean_ub_by_time_CTRL = Df_wholebed_CTRL.groupby('Time')['ub'].mean()
                mean_ub_by_time2 = Df_wholebed2.groupby('Time')['ub'].mean()
                mean_ub_by_time_CTRL2 = Df_wholebed_CTRL2.groupby('Time')['ub'].mean()
                ### Show mean ub in Fig 1
                if LiStyle == '-' and k==0:
                    ax.plot(mean_ub_by_time.index*365, mean_ub_by_time/365, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=LineSize)
                else:
                    ax.plot(mean_ub_by_time.index * 365, mean_ub_by_time / 365, color=blue_col, linestyle=LiStyle, linewidth=LineSize)
                ax.plot(mean_ub_by_time2.index * 365, mean_ub_by_time2 / 365, label= 'C = {}'.format(C), color=blue_col, linestyle='-', linewidth=LineSize)
                ax.plot(mean_ub_by_time_CTRL.index * 365, mean_ub_by_time_CTRL / 365, color=blue_col, linestyle=LiStyle, linewidth=LineSize*0.33)
                ax.plot(mean_ub_by_time_CTRL2.index * 365, mean_ub_by_time_CTRL2 / 365, color=blue_col, linestyle='-', linewidth=LineSize*0.33)
                ###Dummy plots for legend regarding lateral BC
                if i==0:
                    ax.plot(np.nan,np.nan,label='{}'.format(Case),color='k',linestyle=LiStyle, linewidth=LineSize)
                    if j==0 and k==0:
                        ax.plot(np.nan,np.nan,label='CTRL',color='k',linestyle='-', linewidth=LineSize*0.33)
                ax.set_ylim(0.0, 40.0)
                ax.set_xlim(0, 16.01)
                # Major ticks at specific values
                ax.set_xticks([0, 4, 8, 12, 16])##, 20, 24])
                # Minor ticks every 1 unit
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ########################################################
                ###Plot evolution of grounded area over time in Fig2####
                ########################################################
                ###Meaningfull only at bed
                if BoundaryToPlot == 'Bed':
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

                    ####Same thing with new implementation
                    ### Calculate total area over the whole glacier for each time step
                    Total_Area2 = Df_wholebed2.groupby('Time')['nodearea'].sum()
                    Total_Area_CTRL2 = Df_wholebed_CTRL2.groupby('Time')['nodearea'].sum()
                    ### Grounded area
                    Total_GroundedArea2 = (Df_wholebed2[Df_wholebed2['GM'] == 1].groupby('Time')['nodearea'].sum())
                    Total_GroundedArea_CTRL2 = (Df_wholebed_CTRL2[Df_wholebed_CTRL2['GM'] == 1].groupby('Time')['nodearea'].sum())
                    ### Fraction of grounded area over total area
                    Ratio2 = Total_GroundedArea2/Total_Area2
                    Ratio_CTRL2 = Total_GroundedArea_CTRL2/Total_Area_CTRL2
                    ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
                    # Copy ratio but mask values where Ratio == 1
                    Ratio_masked2 = Ratio2.copy()
                    Ratio_masked2[Ratio_masked2 == 1] = np.nan
                    Ratio_masked_CTRL2 = Ratio_CTRL2.copy()
                    Ratio_masked_CTRL2[Ratio_masked_CTRL2 == 1] = np.nan

                    ### Plot ratio
                    if LiStyle == '-':
                        ax2.plot(Total_Area.index*365, Ratio_masked, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=LineSize)
                    else:
                        ax2.plot(Total_Area.index*365, Ratio_masked, color=blue_col,linestyle=LiStyle, linewidth=LineSize)
                    ax2.plot(Total_Area2.index*365, Ratio_masked2, label= 'C = {}'.format(C), color=blue_col,linestyle='-', linewidth=LineSize)
                    ax2.plot(Total_Area_CTRL.index*365, Ratio_masked_CTRL, color=blue_col,linestyle=LiStyle, linewidth=LineSize*0.33)
                    ax2.plot(Total_Area_CTRL2.index*365, Ratio_masked_CTRL2, color=blue_col,linestyle='-', linewidth=LineSize*0.33)
                    ###Dummy plots for legend regarding Case
                    if i==0:
                        ax2.plot(np.nan,np.nan,label='Year = {}'.format(Case),color='k',linestyle=LiStyle, linewidth=LineSize)
                        if j==0 and k==0:
                            ax2.plot(np.nan,np.nan,label='CTRL',color='k',linestyle='-', linewidth=LineSize*0.33)
                    ax2.set_xlim(0, 16.01)
                    # Major ticks at specific values
                    ax2.set_xticks([0, 4, 8, 12, 16])##, 20, 24])
                    # Minor ticks every 1 unit
                    ax2.xaxis.set_minor_locator(MultipleLocator(1))
                    # Add grid for both major and minor ticks
                    ax2.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                    ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

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
    if BoundaryToPlot == 'Bed':
        ax2.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower right')



    plt.show()