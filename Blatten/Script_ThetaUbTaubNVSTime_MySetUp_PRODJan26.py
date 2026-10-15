"""
@author: jbrondex

Description:
------------
This file aims at plotting several variable (taub, theta, theta-dagger, ...) at a prescribed point of the bed
versus time of simu for my setup with progressive loading

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
# Create legend handles manually, for example:
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path
import re
import os
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
import matplotlib.dates as mdates
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib import ticker
from scipy.interpolate import griddata
from matplotlib.ticker import MultipleLocator
import rasterio
from rasterio.windows import from_bounds

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
    # Cvalues = [0.3, 0.4, 0.5] ###Values of C that I want to test
    # LineSizes = [1.8, 2.2, 1.8]
    # LiStyles = [':','-', '--']
    Cvalues = [0.4] ###Values of C that I want to test
    LineSizes = [2.2]
    LiStyles = ['-']
    Year = '2025'
    # StartSimuDay = 16 ####Day of May 2025 corresponding to simulation start
    # SimuDuration = 15 ###Simu duration in days
    StartSimuDay = 23 ####Day of May 2025 corresponding to simulation start
    SimuDuration = 10 ###Simu duration in days
    t0 = pd.Timestamp("2025-05-{} 00:00".format(StartSimuDay))
    tf = t0 + pd.Timedelta(days=SimuDuration)
    CollapseTime = pd.Timestamp("2025-05-28 15:30")
    EndOfLoadingDay = 23
    ###Coordinate of bed point at which we want to plot the variables vs time (this is done in the other python script)
    xplot = [2630580, 2630620.9, 2630653.6]
    yplot = [1139265, 1139221.7, 1139188.3]
    colplot = ['#00FF66','#8D5FD3','#FFD42A']
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 40 ##in m/d !
    Timestep_size = 0.01/365 ### In years !!
    ###Path to Data
    # Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/5_PRONO_PROGRESSIVELOAD_FS_NoGMLat_MaxLoadFromDEMs/GlacierOut/.')
    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/3_PRONO_FULLLOAD_FS_NoGMLat/GlacierOut/.')
    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']

    ###PREPARE THE PLOT
    ###Fig1 is subplot of various variables
    fig1, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    ###Fig 2 same as Fig 1 but organized differently
    fig2, axes2 = plt.subplots(4, 1, figsize=(12, 8), sharex=True)

    ### START OF SCRIPT###
    ### To have a legend for the case considered instead of the upper right panel
    handles=[]
    ###Loop on each C value considered
    for j,(C,LiSize, LiStyle) in enumerate(zip(Cvalues,LineSizes,LiStyles)):
        str_Cvalue = f"{int(C * 10):02d}"
        ###Fill up handles for the common legend regarding C values
        handles.append(Line2D([0], [0], color='black', linestyle=LiStyle, linewidth=LiSize, label='C = {}'.format(str_Cvalue)))
        ####NOW LOAD DATA CORRESPONDING TO PRONOSTIC RUNS WITH OVERLOAD
        print('Dealing wit simu with C =',str_Cvalue, 'and overload finishing on', EndOfLoadingDay, 'of May')
        # Filename = 'MyBirch_RS_C{}_tsp001d_Start{}May_StopLoad{}May_Stop{}May__bed.dat'.format(str_Cvalue, StartSimuDay,EndOfLoadingDay,StartSimuDay + SimuDuration)
        Filename = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp001d_StabFabFSTop_Alpha10m1_FSBotStabStdy__bed.dat'.format(str_Cvalue, StartSimuDay,EndOfLoadingDay,StartSimuDay + SimuDuration)
        ###load file corresponding to simu as dataframe
        print('Opening file:', Filename)
        Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
        Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
        ###Convert absolute time in terms of date
        Df_wholebed['Date'] = t0 + pd.to_timedelta(Df_wholebed['Time'] * 365, unit='D')
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
        ###Deduce ub from u and un
        Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed['Sigma_nn'] = (Df_wholebed['fx'] * Df_wholebed['nx'] + Df_wholebed['fy'] * Df_wholebed['ny'] + Df_wholebed['fz'] * Df_wholebed['nz']) / Df_wholebed['nodearea']
        # List of variables to interpolate
        variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn', 'GM']]
        ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes:
        Df_wholebed.loc[Df_wholebed['GM'] < 0, 'taub'] = 0
        ###Loop on each of the point at which we want a curve
        for i,(xstar,ystar,colstar) in enumerate(zip(xplot,yplot,colplot)):
            # Prepare result container
            results = []
            # Loop over each time step
            for time_val, group in Df_wholebed.groupby('Date'):
                # Extract coordinates of nodes and values
                points = group[['X', 'Y']].values
                # Interpolate each variable at the desired point
                interpolated_values = {}
                for var in variable_columns:
                    values = group[var].values
                    interp_val = griddata(points, values, (xstar, ystar), method='linear')
                    interpolated_values[var] = interp_val
                # Store result
                interpolated_values.update({'Date': time_val, 'X': xstar, 'Y': ystar})
                results.append(interpolated_values)
            # Build final DataFrame
            interpolated_df = pd.DataFrame(results)
            # Sort by time
            interpolated_df = interpolated_df.sort_values(by='Date').reset_index(drop=True)
            interpolated_df .loc[interpolated_df ["GM"] < 0, "taub"] = 0

            #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            #### FROM THERE I START THE PLOT ####
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 1 ####  : Subplots of theta, taub/N, ub, Sigma_nn
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes[0, 0]
            ax.plot(interpolated_df['Date'], (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,linewidth=LiSize)
            ax.plot(interpolated_df['Date'], (1 - interpolated_df['theta_dagger']), color=colstar,linestyle=LiStyle, linewidth=0.8)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\theta$', fontsize=20)
            ax.set_ylim(0.0, 1.0)
            #### DUMMY plot for legend
            if LiStyle == '-' and i == 0:
                ax.plot(np.nan, np.nan, label=r'$\theta$', color='k', linewidth=LiSize, linestyle=LiStyle)
                ax.plot(np.nan, np.nan, label=r'$\theta_{dagger}$', color='k', linewidth=0.8, linestyle=LiStyle)
            ax.legend()
            ax.set_xlim(t0, tf)
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig1.autofmt_xdate()

            ###Subplot 2 = nothing for the moment
            axes[0, 1].axis('off')

            ###Subplot 3 = taub/N
            ax = axes[1, 0]
            ###Plot horizontal line for value of C
            if C <= 0.6:
                ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
                # Add text label above the line on the left side
                ax.text(
                    35, C,  # x and y coordinates in data units
                    'C = {}'.format(C),
                    ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                    fontsize=10, color='gray'
                )
            ax.plot(interpolated_df['Date'], interpolated_df['taub'] / interpolated_df['N'], color=colstar, linestyle=LiStyle, linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.0, 0.6)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)

            ###Subplot 4 = nothing for the moment -> use empty subplot to place legend
            axes[1, 1].axis('off')
            # Place legend in the empty subplot
            axes[1, 1].legend(
                handles=handles,
                loc='center',
                fontsize=14,
                frameon=False  # optional: remove box
            )

            ###Subplot 5 : ub
            ax = axes[2, 0]
            ax.plot(interpolated_df['Date'], interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_ylim(0.0, velocitymax)
            ###Because axes[2,1] is off, put xlabel and force ticks here
            ax.tick_params(labelbottom=True)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)

            ###Subplot 6 : Sigmann
            ax = axes[2, 1]
            ax.plot(interpolated_df['Date'], interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.set_ylim(-5.0, 5.0)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)


            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 2 ####  : Same as Fig 1 but one single column
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes2[0]
            ax.plot(interpolated_df['Date'], (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle, linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\theta$', fontsize=20)
            ax.set_ylim(0.0, 1.0)
            ax.set_xlim(t0, tf)
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)
            fig2.autofmt_xdate()

            ###Subplot 2 = taub/N
            ax = axes2[1]
            ###Plot horizontal line for value of C
            if C <= 0.6:
                ax.axhline(y=C, color='gray', linestyle=LiStyle, linewidth=1.0)
                # Add text label above the line on the left side
                ax.text(
                    16, C,  # x and y coordinates in data units
                    'C = {}'.format(C),
                    ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                    fontsize=10, color='gray'
                )
            ####Everywhere taub is zero I replace by NaN to have a gap in the solid line (makes visualisation easier)
            interpolated_df_taubNaN = interpolated_df.copy()
            interpolated_df_taubNaN.loc[interpolated_df_taubNaN["taub"] == 0, "taub"] = np.nan
            ax.plot(interpolated_df_taubNaN['Date'], interpolated_df_taubNaN['taub'] / interpolated_df_taubNaN['N'], color=colstar,
                    linestyle=LiStyle, linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.0, 0.6)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)


            ###Subplot 3 : ub
            ax = axes2[2]
            ax.plot(interpolated_df['Date'], interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_ylim(0.0, velocitymax)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)

            ###Subplot 4 : Sigmann
            ax = axes2[3]
            ax.plot(interpolated_df['Date'], interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,linewidth=LiSize)
            ###Add vertical line for collapse day
            if i == 0:
                ax.axvline(CollapseTime, linestyle='--', linewidth=2.4, color='r')
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_ylim(-5.0, 5.0)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.6)


    ###Now add legend regarding C on empty upper right subplot
    axes[1, 1].legend(handles=handles,loc='center',fontsize=22,frameon=False)

    plt.show()