"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean ub and grounded area over time with contact problem for all Cs for the geometry of 2025,2011, and 1993

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
# Create legend handles manually, for example:
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.cm as cmx
import re
import os

from pathlib import Path

import numpy as np
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
    Cvalues= [0.3, 0.4, 0.5, 0.6, 0.7]
    # Get 5 shades from the "Blues" colormap (light to dark)
    colors = cmx.get_cmap('viridis_r', len(Cvalues))
    blues = [colors(i) for i in range(len(Cvalues))]
    ###Cases for which we want a figure:
    Cases = ['1993', '2011', '2025']
    LineStyles = [':','--','-']
    # Cases = ['2025']
    # LineStyles = ['-']
    ###Col names
    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
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
        ###Loop on considered years
        for j, (Case, LiStyle) in enumerate(zip(Cases, LineStyles)):
            print('Dealing with year :', Case)
            if Case == '2025':
                Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
                Filename = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue)
                # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
                pattern = re.compile(r"MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(str_Cvalue))
            else:
                Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact_{}/GlacierOut/.'.format(Case))
                Filename = 'MyBirch{}_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d__bed.dat'.format(Case,str_Cvalue)
                # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
                pattern = re.compile(r"MyBirch{}_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(Case,str_Cvalue))
            ##Load data corresponding to first simu
            Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
            ###For year 2025, there is overlap between end of simu R0 and beginning of simu R1
            ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 701 for C=0.3)
            if Case == '2025':
                if C == 0.3:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 701].reset_index(drop=True)
                elif C==0.4:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1101].reset_index(drop=True)
                elif C==0.5:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1601].reset_index(drop=True)
                elif C==0.6:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2501].reset_index(drop=True)
                elif C==0.7:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2901].reset_index(drop=True)
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
            last_timestep = Df_wholebed["Timestep"].max()
            last_count = Df_wholebed["Count"].max()
            ###Now loop on each restart output for current C and year
            for Rstt, RsttPos, fname in matches:
                df_rstt = pd.read_csv(Pathroot.joinpath(fname), names=Col_Name, sep=r'\s+')
                ###For year 2025, there is a problem of overlap between lines of Rj and lines of Rj+1
                if Case == '2025' and Rstt < matches[-1][0]:###This operation is meaningful only if we are not at the last restart
                    ###Get restart position of the next restart
                    NextRstt = Rstt+1
                    ###Get corresponding RsttPos
                    NextRsttPos = next((pos for (rstt, pos, fname) in matches if rstt == NextRstt),
                        None )# default if not found
                    ###Remove overlapping lines
                    df_rstt = df_rstt[df_rstt['Timestep'] <= NextRsttPos+1].reset_index(drop=True)
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
            if Restart_dfs:##Only if Restart_dfs is not empty (there are restarts !)
                Df_wholebed = pd.concat([Df_wholebed] + Restart_dfs, ignore_index=True)
            ###########################################################################
            ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED YEAR####
            Timestep_size = 0.01 / 365  ### In years !!
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
                ax.legend(fontsize=22)
            else:
                ax.plot(mean_ub_by_time.index * 365, mean_ub_by_time / 365, color=blue_col, linestyle=LiStyle, linewidth=2.4)
            ax.set_ylim(0.0, 40.0)
            ax.set_xlim(0, 24.01)
            # Major ticks at specific values
            ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
            # Minor ticks every 1 unit
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ###Now the grounded area for all years
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
                ax2.legend(fontsize=22)
            else:
                ax2.plot(Total_Area.index*365, Ratio_masked, color=blue_col,linestyle=LiStyle, linewidth=2.4)
            ax2.set_xlim(0, 24.01)
            # Major ticks at specific values
            ax2.set_xticks([0, 4, 8, 12, 16, 20, 24])
            # Minor ticks every 1 unit
            ax2.xaxis.set_minor_locator(MultipleLocator(1))
            # Add grid for both major and minor ticks
            ax2.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)




    plt.show()