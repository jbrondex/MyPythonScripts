"""
@author: jbrondex

Description:
------------
This file aims at plotting several variable (taub, theta, theta-dagger, ...) at a prescribed point of the bed
versus time of simu for my setup

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
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
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
    ### In this script I want to compare with and without contact at given C value
    C = 0.3
    str_Cvalue = f"{int(C * 10):02d}"
    str_Year = '2025' ###year of interest FOR NOW WORKS ONLY WITH 2025
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 40 ##in m/d !
    Timestep_size = 0.01/365 ### In years !!
    ###Coordinate of bed point at which we want to plot the variables vs time (this is done in the other python script)
    xplot = [2630580, 2630620.9, 2630653.6]
    yplot = [1139265, 1139221.7, 1139188.3]
    colplot = ['#00FF66','#8D5FD3','#FFD42A']

    ###Path to Data
    Pathroot1 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/2_RateAndState_NoDamage_NoContact/GlacierOut/.')
    Filename1 = 'MyBirch_RS_C{}_WithD_NoCont__bed.dat'.format(str_Cvalue)
    Pathroot2 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
 ###   Filename2 = 'MyBirch_RS_C{}_Cont_SScalInt10_tsp001d_OutVTU1_DVT_R0_FullPicard__bed.dat'.format(str_Cvalue)
    Filename2 = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue)
    Col_Name1 =  ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi']
    Col_Name2 =  ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    ###Line styles for the various case
    Pathroots= [Pathroot2]#[Pathroot1, Pathroot2]
    Filenames = [Filename2]#[Filename1, Filename2]
    Col_Names = [Col_Name2]#[Col_Name1, Col_Name2]
    LiStyles = ['-']#['--', '-']
    Labels = ['Contact']#['No Contact', 'Contact']

    ###PREPARE THE PLOT
    ###Fig1 is subplot of various variables
    fig1, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    ###Fig 2 same as Fig 1 but organized differently
    fig2, axes2 = plt.subplots(4, 1, figsize=(12, 8), sharex=True)
    ####Figure 3 : tb/N=f(ub/N)
    # Create dedicated figure and axis
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    ###Fig 4 same as Fig 2 but comparing FreeSurface Lateral versus no normal velo for C =0.3 only
    if C ==0.3:
        fig4, axes4 = plt.subplots(4, 1, figsize=(12, 8), sharex=True)

    ### START OF SCRIPT###
    ### To have a legend for the case considered instead of the upper right panel
    handles=[]
    ###Loop on each case (contact vs no contact) considered
    for j,(Pathroot, Filename, Col_Name , Label, LiStyle) in enumerate(zip(Pathroots, Filenames, Col_Names, Labels, LiStyles)):
        ###Correct time step size for old simu where I was using tsp 0.02 days:
        if Filename == 'MyBirch_RS_C{}_WithD_NoCont__bed.dat'.format(str_Cvalue):
            Timestep_size = 0.02/365 ### In years !!
        ###Fill up handles for the common legend regarding C values
        handles.append(Line2D([0], [0], color='black', linestyle=LiStyle, linewidth=3, label=Label))
        ###load file as dataframe
        Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
        ####For the case with contact problem solved, there are restarts to handle
        if Filename ==  'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue):
            # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
            pattern = re.compile(r"MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(str_Cvalue))
            if str_Year== '2025': ###COMPULSORY FOR NOW
                ###For year 2025, there is overlap between end of simu R0 and beginning of simu R1
                ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 701 for C=0.3)
                if C == 0.3:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 701].reset_index(drop=True)
                elif C == 0.4:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1101].reset_index(drop=True)
                elif C == 0.5:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1601].reset_index(drop=True)
                elif C == 0.6:
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2501].reset_index(drop=True)
                elif C == 0.7:
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
                if Rstt < matches[-1][0]:  ###This operation is meaningful only if we are not at the last restart
                    ###Get restart position of the next restart
                    NextRstt = Rstt + 1
                    ###Get corresponding RsttPos
                    NextRsttPos = next((pos for (rstt, pos, fname) in matches if rstt == NextRstt),None)  # default if not found
                    ###Remove overlapping lines
                    df_rstt = df_rstt[df_rstt['Timestep'] <= NextRsttPos + 1].reset_index(drop=True)
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
                Df_wholebed = pd.concat([Df_wholebed] + Restart_dfs, ignore_index=True)

        ###Create a column with time (in years !)
        Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
        # We look only the 40 first days : Keep only rows where Time <= 40/365 and reset index
        Df_wholebed = Df_wholebed[Df_wholebed['Time'] <= 40 / 365].reset_index(drop=True)
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
        ###Deduce ub from u and un
        Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed['Sigma_nn'] = (Df_wholebed['fx'] * Df_wholebed['nx'] + Df_wholebed['fy'] * Df_wholebed['ny'] + Df_wholebed['fz'] * Df_wholebed['nz']) / Df_wholebed['nodearea']
        # List of variables to interpolate
        # variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'Chi', 'ub', 'Sigma_nn']]
        if Pathroot == Pathroot1: ##For the no contact case, no GM
            variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn']]
        elif Pathroot == Pathroot2: ##For the with contact case, we also interpolate the GM to correct the diagnosed taub
            variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn', 'GM']]
            ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes:
            Df_wholebed.loc[Df_wholebed['GM'] < 0, 'taub'] = 0
        ###Loop on each of the point at which we want a curve
        for i,(xstar,ystar,colstar) in enumerate(zip(xplot,yplot,colplot)):
            # Prepare result container
            results = []
            # Loop over each time step
            for time_val, group in Df_wholebed.groupby('Time'):
                # Extract coordinates of nodes and values
                points = group[['X', 'Y']].values
                # Interpolate each variable at the desired point
                interpolated_values = {}
                for var in variable_columns:
                    values = group[var].values
                    interp_val = griddata(points, values, (xstar, ystar), method='linear')
                    interpolated_values[var] = interp_val
                # Store result
                interpolated_values.update({'Time': time_val, 'X': xstar, 'Y': ystar})
                results.append(interpolated_values)
            # Build final DataFrame
            interpolated_df = pd.DataFrame(results)
            # Sort by time
            interpolated_df = interpolated_df.sort_values(by='Time').reset_index(drop=True)
            # # ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes for the with contact case:
            if Pathroot == Pathroot2: ##Only for the case with contact
                # Force taub = 0 where GM < 0
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
            ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,
                    linewidth=2.2)
            ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta_dagger']), color=colstar,
                    linestyle=LiStyle, linewidth=0.8)
            ax.set_ylabel(r'$\theta$', fontsize=20)
            ax.set_ylim(0.0, 1.0)
            #### DUMMY plot for legend
            if LiStyle == '-' and i == 0:
                ax.plot(np.nan, np.nan, label=r'$\theta$', color='k', linewidth=2.2, linestyle=LiStyle)
                ax.plot(np.nan, np.nan, label=r'$\theta_{dagger}$', color='k', linewidth=0.8, linestyle=LiStyle)
            ax.legend()
            ax.grid(True)

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
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['taub'] / interpolated_df['N'], color=colstar,
                    linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.2, 0.6)
            ax.grid(True)

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
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_ylim(0.0, velocitymax)
            ax.grid(True)
            ###Because axes[2,1] is off, put xlabel and force ticks here
            ax.tick_params(labelbottom=True)
            ax.set_xlabel(r'Time [days]', fontsize=20)

            ###Subplot 6 : Sigmann
            ax = axes[2, 1]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.set_ylim(-5.0, 5.0)
            ax.grid(True)

            ####Spare Subplot : Chi
            # ax.plot(interpolated_df['Time']*365, interpolated_df['Chi'], color=colstar, linestyle=LiStyle, linewidth=2)
            # ax.axhspan(0, 1.8, color='grey', alpha=0.05)
            # ax.set_ylabel(r'$\chi$', fontsize=20)
            # ax.set_xlabel(r'Time [days]', fontsize=20)
            # ax.set_ylim(-3.5,1.8)
            # ax.grid(True)

            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 2 ####  : Same as Fig 1 but one single column
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes2[0]
            ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,
                    linewidth=2.2)
            ax.set_ylabel(r'$\theta$', fontsize=20)
            ax.set_ylim(0.0, 1.0)
            ###xlim for all plots : only first 20 days
            ax.set_xlim(0, 24.01)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ###Subplot 2 = taub/N
            ax = axes2[1]
            ###Plot horizontal line for value of C
            if C <= 0.6:
                ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
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
            ax.plot(interpolated_df_taubNaN['Time'] * 365, interpolated_df_taubNaN['taub'] / interpolated_df_taubNaN['N'], color=colstar,
                    linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.2, 0.7)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)


            ###Subplot 3 : ub
            ax = axes2[2]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_ylim(0.0, velocitymax)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ###Subplot 4 : Sigmann
            ax = axes2[3]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.set_ylim(-5.0, 5.0)
            # Major ticks at specific values
            ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
            # Minor ticks every 1 unit
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 3 ####
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            # Fill up figure 3 : Plot the ratio ub/N vs taub/N
            ax3.plot(interpolated_df['ub'] / (365 * interpolated_df['N']),
                     interpolated_df['taub'] / interpolated_df['N'], color=colstar, linestyle=LiStyle, linewidth=2)
            # Axis labels with LaTeX-style formatting
            ax3.set_xlabel(r'$\mathrm{u}_\mathrm{b} / N$ [m d$^{-1}$ MPa$^{-1}$]', fontsize=20)
            ax3.set_ylabel(r'$\tau_b / N$', fontsize=20)
            ax3.grid(True)

            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 4 ####  : Same as Fig 2 but for the case of free surf lat (C=0.3 only)
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            if C == 0.3:
                ###Subplot 1 : Theta
                ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
                ax = axes4[0]
                ###xlim for all plots : only first 20 days
                ax.set_xlim(0, 24.01)
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,
                        linewidth=2.2)
                ax.set_ylabel(r'$\theta$', fontsize=20)
                ax.set_ylim(0.0, 1.0)



                ###Subplot 2 = taub/N
                ax = axes4[1]
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ###Plot horizontal line for value of C
                if C <= 0.6:
                    ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
                    # Add text label above the line on the left side
                    ax.text(
                        21, C,  # x and y coordinates in data units
                        'C = {}'.format(C),
                        ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                        fontsize=14, color='gray'
                    )
                ####Everywhere taub is zero I replace by NaN to have a gap in the solid line (makes visualisation easier)
                interpolated_df_taubNaN = interpolated_df.copy()
                interpolated_df_taubNaN.loc[interpolated_df_taubNaN["taub"] == 0, "taub"] = np.nan
                ax.plot(interpolated_df_taubNaN['Time'] * 365, interpolated_df_taubNaN['taub'] / interpolated_df_taubNaN['N'], color=colstar,
                        linestyle=LiStyle, linewidth=2)
                ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
                ax.set_ylim(0.2, 1.0)


                ###Subplot 3 : ub
                ax = axes4[2]
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ax.plot(interpolated_df['Time'] * 365, interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,
                        linewidth=2)
                ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
                ax.set_ylim(0.0, velocitymax+10)


                ###Subplot 4 : Sigmann
                ax = axes4[3]
                # Major ticks at specific values
                ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
                # Minor ticks every 1 unit
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ax.plot(interpolated_df['Time'] * 365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,
                        linewidth=2)
                ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
                ax.set_xlabel(r'Time [days]', fontsize=20)
                ax.set_ylim(-7.0, 2.2)


    #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    #### HERE I TAKE CARE OF THE CASE FREE SURF LATERAL FOR C =0.3 and WITH CONTACT PROBLEM ONLY
    if C == 0.3:
        Filename3 = 'MyBirch_RS_C03_Cont_SScalInt10_tsp001d_OutVTU1_DVT_FullPicard_LatFS__bed.dat'
        Timestep_size = 0.01 / 365  ### In years !!
        ###Fill up handles for the common legend regarding C values
        Df_wholebed3 = pd.DataFrame()
        ###load file as dataframe
        Df_wholebed3 = pd.read_csv(Pathroot2.joinpath(Filename3), names=Col_Name2, sep=r'\s+')
        ###Create a column with time (in years !)
        Df_wholebed3['Time'] = Timestep_size * Df_wholebed3['Timestep']
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed3['un'] = Df_wholebed3['u'] * Df_wholebed3['nx'] + Df_wholebed3['v'] * Df_wholebed3['ny'] + \
                             Df_wholebed3['w'] * Df_wholebed3['nz']
        ###Deduce ub from u and un
        Df_wholebed3['ub'] = np.sqrt((Df_wholebed3['u'] - Df_wholebed3['un'] * Df_wholebed3['nx']) ** 2 + (
                Df_wholebed3['v'] - Df_wholebed3['un'] * Df_wholebed3['ny']) ** 2 + (
                                             Df_wholebed3['w'] - Df_wholebed3['un'] * Df_wholebed3['nz']) ** 2)
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed3['Sigma_nn'] = (Df_wholebed3['fx'] * Df_wholebed3['nx'] + Df_wholebed3['fy'] * Df_wholebed3[
            'ny'] + Df_wholebed3['fz'] * Df_wholebed3['nz']) / Df_wholebed3['nodearea']
        # List of variables to interpolate
        variable_columns = [col for col in Df_wholebed3.columns if
                            col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn', 'GM']]
        ###Loop on each of the point at which we want a curve
        for i, (xstar, ystar, colstar) in enumerate(zip(xplot, yplot, colplot)):
            # Prepare result container
            results3 = []
            # Loop over each time step
            for time_val, group in Df_wholebed3.groupby('Time'):
                # Extract coordinates of nodes and values
                points = group[['X', 'Y']].values
                # Interpolate each variable at the desired point
                interpolated_values3 = {}
                for var in variable_columns:
                    values = group[var].values
                    interp_val = griddata(points, values, (xstar, ystar), method='linear')
                    interpolated_values3[var] = interp_val
                # Store result
                interpolated_values3.update({'Time': time_val, 'X': xstar, 'Y': ystar})
                results3.append(interpolated_values3)
            # Build final DataFrame
            interpolated_df3 = pd.DataFrame(results3)
            # Sort by time
            interpolated_df3 = interpolated_df3.sort_values(by='Time').reset_index(drop=True)
            # Force taub = 0 where GM < 0
            interpolated_df3 .loc[interpolated_df3 ["GM"] < 0, "taub"] = 0
            ### FILL-UP FIGURE 4 WITH THIS LATERAL FREE SURFACE CASE
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes4[0]
            ax.plot(interpolated_df3['Time'] * 365, (1 - interpolated_df3['theta']), color=colstar, linestyle=':',linewidth=2.2)
            ###Subplot 2 = taub/N
            ax = axes4[1]
            ####Everywhere taub is zero I replace by NaN to have a gap in the solid line (makes visualisation easier)
            interpolated_df3_taubNaN = interpolated_df3.copy()
            interpolated_df3_taubNaN.loc[interpolated_df3_taubNaN["taub"] == 0, "taub"] = np.nan
            ax.plot(interpolated_df3_taubNaN['Time'] * 365, interpolated_df3_taubNaN['taub'] / interpolated_df3_taubNaN['N'], color=colstar,linestyle=':', linewidth=2)
            ###Subplot 3 : ub
            ax = axes4[2]
            ax.plot(interpolated_df3['Time'] * 365, interpolated_df3['ub'] / 365, color=colstar, linestyle=':',linewidth=2)
            ###Subplot 4 : Sigmann
            ax = axes4[3]
            ax.plot(interpolated_df3['Time'] * 365, interpolated_df3['Sigma_nn'], color=colstar, linestyle=':',linewidth=2)

    ###Now add legend regarding C on empty upper right subplot
    axes[1, 1].legend(handles=handles,loc='center',fontsize=22,frameon=False)

    plt.show()