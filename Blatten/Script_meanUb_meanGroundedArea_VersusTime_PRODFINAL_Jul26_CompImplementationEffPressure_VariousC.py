"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean ub and grounded area over time for chosen C values to compare results for two implementations
of the rate and state friction law : (i) the older one where effective pressure is computed assuming that normal stress of ice on bed
is given by ice overburden pressure, (ii) the most recent one (July 2026) where effective pressure is computed from the true normal stress at bed

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
    Cvalues= [0.3, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.4]
    # Get 5 shades from the "Blues" colormap (light to dark)
    colors = cmx.get_cmap('viridis_r', len(Cvalues))
    blues = [colors(i) for i in range(len(Cvalues))]
    LineSize = 2.4 ###Default value of linewidth

    Pathroot_May26 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
    Pathroot_July26 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')

    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    Col_Name2 = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'Normal Stress', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax.set_xlabel(r'Time [days]', fontsize=20)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=20)
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
        ####Rearrangement of parameters to get corresponding file names
        s = f"{C:.2f}"
        if len(s.rstrip('0').split('.')[1]) == 1:
            str_Cvalue = f"{int(C * 10):02d}"
        elif len(s.rstrip('0').split('.')[1]) == 2:
            str_Cvalue = f"{int(C * 100):03d}"
        Filename_May26 = 'MyBirch_RS_C{}_tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(str_Cvalue)
        Filename_July26 = 'MyBirch_RS_C{}_dc2m_tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(str_Cvalue)
        ###Time step size for prono simu with overload:
        Timestep_size = 0.01 / 365  ### In years !!
        ###load file corresponding to first simu as dataframe
        print('Opening file:', Filename_May26)
        Df_wholebed_May26 = pd.read_csv(Pathroot_May26.joinpath(Filename_May26), names=Col_Name, sep=r'\s+')
        Df_wholebed_May26 = Df_wholebed_May26[Df_wholebed_May26['BC'] == 3]
        print('Opening file:', Filename_July26)
        Df_wholebed_July26 = pd.read_csv(Pathroot_July26.joinpath(Filename_July26), names=Col_Name2, sep=r'\s+')
        Df_wholebed_July26 = Df_wholebed_July26[Df_wholebed_July26['BC'] == 3]

        ###########################################################################
        ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED CASE####
        ###Create a column with time (in years !)
        Df_wholebed_May26['Time'] = Timestep_size * Df_wholebed_May26['Timestep']
        Df_wholebed_July26['Time'] = Timestep_size * Df_wholebed_July26['Timestep']

        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed_May26['un'] = Df_wholebed_May26['u'] * Df_wholebed_May26['nx'] + Df_wholebed_May26['v'] * Df_wholebed_May26['ny'] + Df_wholebed_May26['w'] * Df_wholebed_May26['nz']
        ###Deduce ub from u and un
        Df_wholebed_May26['ub'] = np.sqrt((Df_wholebed_May26['u'] - Df_wholebed_May26['un'] * Df_wholebed_May26['nx']) ** 2 + (Df_wholebed_May26['v'] - Df_wholebed_May26['un'] * Df_wholebed_May26['ny']) ** 2 + (Df_wholebed_May26['w'] - Df_wholebed_May26['un'] * Df_wholebed_May26['nz']) ** 2)
        ####Samething with new implementation
        Df_wholebed_July26 ['un'] = Df_wholebed_July26 ['u'] * Df_wholebed_July26 ['nx'] + Df_wholebed_July26 ['v'] * Df_wholebed_July26 ['ny'] + Df_wholebed_July26 ['w'] * Df_wholebed_July26 ['nz']
        ###Deduce ub from u and un
        Df_wholebed_July26 ['ub'] = np.sqrt((Df_wholebed_July26 ['u'] - Df_wholebed_July26 ['un'] * Df_wholebed_July26 ['nx']) ** 2 + (Df_wholebed_July26 ['v'] - Df_wholebed_July26 ['un'] * Df_wholebed_July26 ['ny']) ** 2 + (Df_wholebed_July26 ['w'] - Df_wholebed_July26 ['un'] * Df_wholebed_July26 ['nz']) ** 2)

        ### Calculate mean ub over whole glacier for each time step
        mean_ub_by_time_May26 = Df_wholebed_May26.groupby('Time')['ub'].mean()
        mean_ub_by_time_July26 = Df_wholebed_July26.groupby('Time')['ub'].mean()
        ###DO THE PLOT
        ax.plot(mean_ub_by_time_May26.index*365, mean_ub_by_time_May26/365, color=blue_col,linestyle='-', linewidth=LineSize*0.33)
        ax.plot(mean_ub_by_time_July26.index*365, mean_ub_by_time_July26/365, label= 'C = {}'.format(C), color=blue_col,linestyle='-', linewidth=LineSize)
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

        ### Calculate total area over the whole glacier for each time step
        Total_Area_May26 = Df_wholebed_May26.groupby('Time')['nodearea'].sum()
        ### Grounded area
        Total_GroundedArea_May26 = (Df_wholebed_May26[Df_wholebed_May26['GM'] == 1].groupby('Time')['nodearea'].sum())
        ### Fraction of grounded area over total area
        Ratio_May26 = Total_GroundedArea_May26/Total_Area_May26
        ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
        # Copy ratio but mask values where Ratio == 1
        Ratio_masked_May26 = Ratio_May26.copy()
        Ratio_masked_May26[Ratio_masked_May26 == 1] = np.nan

        ####Same thing with new implementation
        ### Calculate total area over the whole glacier for each time step
        Total_Area_July26 = Df_wholebed_July26.groupby('Time')['nodearea'].sum()
        ### Grounded area
        Total_GroundedArea_July26 = (Df_wholebed_July26[Df_wholebed_July26['GM'] == 1].groupby('Time')['nodearea'].sum())
        ### Fraction of grounded area over total area
        Ratio_July26 = Total_GroundedArea_July26/Total_Area_July26
        ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
        # Copy ratio but mask values where Ratio == 1
        Ratio_masked_July26 = Ratio_July26.copy()
        Ratio_masked_July26[Ratio_masked_July26 == 1] = np.nan

        ### Plot ratio
        ax2.plot(Total_Area_May26.index*365, Ratio_masked_May26, color=blue_col,linestyle='-', linewidth=LineSize*0.33)
        ax2.plot(Total_Area_July26.index*365, Ratio_masked_July26, label= 'C = {}'.format(C), color=blue_col,linestyle='-', linewidth=LineSize)
        ###Dummy plots for legend regarding Case
        ax2.set_xlim(0, 16.01)
        # Major ticks at specific values
        ax2.set_xticks([0, 4, 8, 12, 16])##, 20, 24])
        # Minor ticks every 1 unit
        ax2.xaxis.set_minor_locator(MultipleLocator(1))
        # Add grid for both major and minor ticks
        ax2.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
        ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

    ### NOW ADD CASES FOR WHICH DC IS NOT 2M ####
    ####For C=0.36 we want to plot also the case with dc=1.5m
    ####For C=0.3 and C=0.35 we want to plot the case dc=3.0m
    for j, (C, dc, LiStyle) in enumerate(zip([0.36], [1.5], ['--'])):
        colorplot = blues[Cvalues.index(C)]
        ####Rearrangement of parameters to get corresponding file names
        s = f"{C:.2f}"
        if len(s.rstrip('0').split('.')[1]) == 1:
            str_Cvalue = f"{int(C * 10):02d}"
        elif len(s.rstrip('0').split('.')[1]) == 2:
            str_Cvalue = f"{int(C * 100):03d}"

        if dc.is_integer():
            str_dcvalue = f"dc{int(dc):01d}m_"
        else:
            x, y = str(dc).split(".")
            str_dcvalue = f"dc{x}dot{y}m_"
        Filename_May26 = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue)
        Filename_July26 = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue)
        ###Time step size for prono simu with overload:
        Timestep_size = 0.01 / 365  ### In years !!
        ###load file corresponding to first simu as dataframe
        print('Opening file:', Filename_May26)
        Df_wholebed_May26 = pd.read_csv(Pathroot_May26.joinpath(Filename_May26), names=Col_Name, sep=r'\s+')
        Df_wholebed_May26 = Df_wholebed_May26[Df_wholebed_May26['BC'] == 3]
        print('Opening file:', Filename_July26)
        Df_wholebed_July26 = pd.read_csv(Pathroot_July26.joinpath(Filename_July26), names=Col_Name2, sep=r'\s+')
        Df_wholebed_July26 = Df_wholebed_July26[Df_wholebed_July26['BC'] == 3]

        ###########################################################################
        ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED CASE####
        ###Create a column with time (in years !)
        Df_wholebed_May26['Time'] = Timestep_size * Df_wholebed_May26['Timestep']
        Df_wholebed_July26['Time'] = Timestep_size * Df_wholebed_July26['Timestep']

        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed_May26['un'] = Df_wholebed_May26['u'] * Df_wholebed_May26['nx'] + Df_wholebed_May26['v'] * \
                                  Df_wholebed_May26['ny'] + Df_wholebed_May26['w'] * Df_wholebed_May26['nz']
        ###Deduce ub from u and un
        Df_wholebed_May26['ub'] = np.sqrt(
            (Df_wholebed_May26['u'] - Df_wholebed_May26['un'] * Df_wholebed_May26['nx']) ** 2 + (
                        Df_wholebed_May26['v'] - Df_wholebed_May26['un'] * Df_wholebed_May26['ny']) ** 2 + (
                        Df_wholebed_May26['w'] - Df_wholebed_May26['un'] * Df_wholebed_May26['nz']) ** 2)
        ####Samething with new implementation
        Df_wholebed_July26['un'] = Df_wholebed_July26['u'] * Df_wholebed_July26['nx'] + Df_wholebed_July26['v'] * \
                                   Df_wholebed_July26['ny'] + Df_wholebed_July26['w'] * Df_wholebed_July26['nz']
        ###Deduce ub from u and un
        Df_wholebed_July26['ub'] = np.sqrt(
            (Df_wholebed_July26['u'] - Df_wholebed_July26['un'] * Df_wholebed_July26['nx']) ** 2 + (
                        Df_wholebed_July26['v'] - Df_wholebed_July26['un'] * Df_wholebed_July26['ny']) ** 2 + (
                        Df_wholebed_July26['w'] - Df_wholebed_July26['un'] * Df_wholebed_July26['nz']) ** 2)

        ### Calculate mean ub over whole glacier for each time step
        mean_ub_by_time_May26 = Df_wholebed_May26.groupby('Time')['ub'].mean()
        mean_ub_by_time_July26 = Df_wholebed_July26.groupby('Time')['ub'].mean()
        ###DO THE PLOT
        ax.plot(mean_ub_by_time_May26.index * 365, mean_ub_by_time_May26 / 365, color=colorplot, linestyle=LiStyle,
                linewidth=LineSize * 0.33)
        ax.plot(mean_ub_by_time_July26.index * 365, mean_ub_by_time_July26 / 365, color=colorplot,
                linestyle=LiStyle, linewidth=LineSize)
        ########################################################
        ###Plot evolution of grounded area over time in Fig2####
        ########################################################
        ### Calculate total area over the whole glacier for each time step
        Total_Area_May26 = Df_wholebed_May26.groupby('Time')['nodearea'].sum()
        ### Grounded area
        Total_GroundedArea_May26 = (Df_wholebed_May26[Df_wholebed_May26['GM'] == 1].groupby('Time')['nodearea'].sum())
        ### Fraction of grounded area over total area
        Ratio_May26 = Total_GroundedArea_May26 / Total_Area_May26
        ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
        # Copy ratio but mask values where Ratio == 1
        Ratio_masked_May26 = Ratio_May26.copy()
        Ratio_masked_May26[Ratio_masked_May26 == 1] = np.nan

        ####Same thing with new implementation
        ### Calculate total area over the whole glacier for each time step
        Total_Area_July26 = Df_wholebed_July26.groupby('Time')['nodearea'].sum()
        ### Grounded area
        Total_GroundedArea_July26 = (Df_wholebed_July26[Df_wholebed_July26['GM'] == 1].groupby('Time')['nodearea'].sum())
        ### Fraction of grounded area over total area
        Ratio_July26 = Total_GroundedArea_July26 / Total_Area_July26
        ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
        # Copy ratio but mask values where Ratio == 1
        Ratio_masked_July26 = Ratio_July26.copy()
        Ratio_masked_July26[Ratio_masked_July26 == 1] = np.nan

        ### Plot ratio
        ax2.plot(Total_Area_May26.index * 365, Ratio_masked_May26, color=colorplot, linestyle=LiStyle, linewidth=LineSize * 0.33)
        ax2.plot(Total_Area_July26.index * 365, Ratio_masked_July26, color=colorplot,linestyle=LiStyle, linewidth=LineSize)


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
    ax2.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower left')



    plt.show()