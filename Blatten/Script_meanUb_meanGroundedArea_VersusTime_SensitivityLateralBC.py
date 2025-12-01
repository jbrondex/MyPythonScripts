"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean ub and grounded area over time for all some C with different BCs on the lateral Boundary

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
    Cvalues= [0.3, 0.5]
    # Get 5 shades from the "Blues" colormap (light to dark)
    colors = cmx.get_cmap('viridis_r', len(Cvalues))
    blues = [colors(i) for i in range(len(Cvalues))]
    ###Cases in terms of Lateral BC for which we want a figure:
    Cases = ['GroundedMaskLateralUpdated','GroundedMaskLateral','FreeSurf', 'FreeSlipNoNormalV']#, 'NoVelo']
    LineStyles = ['-','--',':','-.']#,':']
    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub for both with and without contact problem ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=20)
    ax.set_xlabel(r'Time [days]', fontsize=20)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ax.set_title('Mean basal velocity [m/day]', fontsize=21, weight='bold')
    ###~~~~~~~FIG 2 : grounded fraction for with contact problem only ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig2, ax2 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax2.set_ylabel(r'Grounded Area', fontsize=20)
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
        ###Loop on lateral BC considered
        for j, (Case, LiStyle) in enumerate(zip(Cases, LineStyles)):
            print('Dealing with Lateral BC :', Case)
            Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
            Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
            Filename = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_LatBC{}__bed.dat'.format(str_Cvalue,Case)
            # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
            pattern = re.compile(r"MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d__LatBC{}_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(str_Cvalue,Case))
            ###Correct time step size for the simu for which I changes to 0.01 day:
            Timestep_size = 0.01 / 365  ### In years !!
            ###load file as dataframe
            Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
            ###For the case with contact, restarts are required : make the full dataframe from all restarts
            # # Collect restart DataFrames into a list
            # Restart_dfs = []
            # # Gather all restart files and sort them by Rstt (restart number)
            # matches = []
            # for fname in os.listdir(Pathroot):
            #     match = pattern.match(fname)
            #     if match:
            #         Rstt, RsttPos = map(int, match.groups())
            #         matches.append((Rstt, RsttPos, fname))
            # matches.sort(key=lambda x: x[0])  # ensure correct order : Restart 1, then restart 2 and so on
            # # Loop efficiently over sorted restart files
            # last_timestep = Df_wholebed["Timestep"].max()
            # last_count = Df_wholebed["Count"].max()
            # ###Now loop on each restart output for current C and year
            # for Rstt, RsttPos, fname in matches:
            #     df_rstt = pd.read_csv(Pathroot.joinpath(fname), names=Col_Name, sep=r'\s+')
            #     ###For year 2025, there is a problem of overlap between lines of Rj and lines of Rj+1
            #     if Rstt < matches[-1][0]:  ###This operation is meaningful only if we are not at the last restart
            #         ###Get restart position of the next restart
            #         NextRstt = Rstt + 1
            #         ###Get corresponding RsttPos
            #         NextRsttPos = next((pos for (rstt, pos, fname) in matches if rstt == NextRstt),None)  # default if not found
            #         ###Remove overlapping lines
            #         df_rstt = df_rstt[df_rstt['Timestep'] <= NextRsttPos + 1].reset_index(drop=True)
            #     ###NOW TREAT THE DATAFRAME OF THE RESTART
            #     # Skip first timestep and shift
            #     df_rstt = df_rstt[df_rstt["Timestep"] != 1].copy()
            #     df_rstt["Timestep"] += last_timestep - 1
            #     df_rstt["Count"] += last_count - 1
            #     # Update counters
            #     last_timestep = df_rstt["Timestep"].max()
            #     last_count = df_rstt["Count"].max()
            #     ## Add to the full dataframe (first as a list)
            #     Restart_dfs.append(df_rstt)

            # # Concatenate once at the end (much faster)
            # if Restart_dfs:  ##Only if Restart_dfs is not empty (there are restarts !)
            #     Df_wholebed = pd.concat([Df_wholebed] + Restart_dfs, ignore_index=True)
                ###########################################################################
                ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED BC####
            ###Create a column with time (in years !)
            Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
            ###Calculate ub from velocity and normal vector
            ###Start by calculating un (should be zero but there is some residual)
            Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
            ###Deduce ub from u and un
            Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
            ### Calculate mean ub over whole glacier for each time step
            mean_ub_by_time = Df_wholebed.groupby('Time')['ub'].mean()
            ### Show mean ub in Fig 1
            if LiStyle == '-':
                ax.plot(mean_ub_by_time.index*365, mean_ub_by_time/365, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=2.4)
            else:
                ax.plot(mean_ub_by_time.index * 365, mean_ub_by_time / 365, color=blue_col, linestyle=LiStyle, linewidth=2.4)
            ###Dummy plots for legend regarding lateral BC
            if i==0:
                ax.plot(np.nan,np.nan,label='BC = {}'.format(Case),color='k',linestyle=LiStyle, linewidth=2.4)
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
            ### Grounded area
            Total_GroundedArea = (Df_wholebed[Df_wholebed['GM'] == 1].groupby('Time')['nodearea'].sum())
            ### Fraction of grounded area over total area
            Ratio = Total_GroundedArea/Total_Area
            ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
            # Copy ratio but mask values where Ratio == 1
            Ratio_masked = Ratio.copy()
            Ratio_masked[Ratio_masked == 1] = np.nan

            ### Plot ratio
            if LiStyle == '-':
                ax2.plot(Total_Area.index*365, Ratio_masked, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=2.4)
            else:
                ax2.plot(Total_Area.index*365, Ratio_masked, color=blue_col,linestyle=LiStyle, linewidth=2.4)
            ###Dummy plots for legend regarding lateral BC
            if i==0:
                ax2.plot(np.nan,np.nan,label='BC = {}'.format(Case),color='k',linestyle=LiStyle, linewidth=2.4)
            ax2.set_xlim(0, 24.01)
            # Major ticks at specific values
            ax2.set_xticks([0, 4, 8, 12, 16, 20, 24])
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

    for h, l in zip(handles, labels):
        if l.startswith("C ="):
            handles_C.append(h)
            labels_C.append(l)
        elif l.startswith("BC ="):
            handles_BC.append(h)
            labels_BC.append(l)

    # Combine them: C first, BC after
    handles_ordered = handles_C + handles_BC
    labels_ordered = labels_C + labels_BC

    ax.legend(handles_ordered, labels_ordered, fontsize=22)
    ax2.legend(handles_ordered, labels_ordered, fontsize=22)



    plt.show()