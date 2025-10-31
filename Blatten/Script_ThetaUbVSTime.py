"""
@author: jbrondex

Description:
------------
This file aims at plotting several variable (taub, theta, theta-dagger, ...) at a prescribed point of the bed
versus time of simu.

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
    Case = "WithContact" ####Do we want to plot case "WithContact" or "NoConact" ?
    ###Open Data corresponding to grounded mask
    if Case =="WithContact":
        Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
        Col_Names = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    elif Case == "NoContact":
        Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/2_RateAndState_NoDamage_NoContact/GlacierOut/.')
        Col_Names = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi']
    else:
        print(' Variable "Case" at beginning of script should be set to "WithContact" or "NoContact"')
    ###Parameter of Simu
    Timestep_size = 0.02/365 ### In years !!
    Cvalues= [0.3, 0.4, 0.5]#, 0.7]
    LiStyles = [':','--','-']#,'-.']
    ###Coordinate of bed point at which we want to plot the variables vs time (this is done in the other python script)
    xplot = [2630580, 2630620.9, 2630653.6]
    yplot = [1139265, 1139221.7, 1139188.3]
    colplot = ['#00FF66','#8D5FD3','#FFD42A']

    ###PREPARE THE PLOT
    ###Fig1 is subplot of various variables
    fig1, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    ####Figure 2 : tb/N=f(ub/N)
    # Create dedicated figure and axis
    fig2, ax2 = plt.subplots(figsize=(8, 6))

    ### START OF SCRIPT###
    ### To have a legend for the C considered instead of the upper right panel
    handles=[]
    ###Loop on each C considered
    for j,(C,LiStyle) in enumerate(zip(Cvalues, LiStyles)):
        str_Cvalue = f"{int(C * 10):02d}"
        # Filename = 'Birch_rock_rate_C{}_WithDamage_MPS_ExportChi__bed.dat'.format(str_Cvalue)
        if Case == "WithContact":
            Filename = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue)
        elif Case == "NoContact":
            Filename = 'MyBirch_RS_C{}_WithD_NoCont__bed.dat'.format(str_Cvalue)
        else:
            print(' Variable "Case" at beginning of script should be set to "WithContact" or "NoContact"')
        ###Fill up handles for the common legend regarding C values
        handles.append(Line2D([0], [0], color='black', linestyle=LiStyle, linewidth=3, label=f"C = {C}"))
        ###load file as dataframe
        Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Names, sep=r'\s+')
        ###Create a column with time (in years !)
        Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
        ###Deduce ub from u and un
        Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed['Sigma_nn'] = (Df_wholebed['fx'] * Df_wholebed['nx'] + Df_wholebed['fy'] * Df_wholebed['ny'] + Df_wholebed['fz'] * Df_wholebed['nz']) / Df_wholebed['nodearea']
        # List of variables to interpolate
        variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'Chi', 'ub', 'Sigma_nn']]
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
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes[0,0]
            ax.plot(interpolated_df['Time']*365, (1-interpolated_df['theta']), color=colstar, linestyle=LiStyle, linewidth=2.2)
            ax.plot(interpolated_df['Time']*365, (1-interpolated_df['theta_dagger']), color=colstar, linestyle=LiStyle, linewidth=0.8)
            ax.set_ylabel(r'$\theta$', fontsize=20)
            #### DUMMY plot for legend
            if LiStyle=='-' and i ==0:
                ax.plot(np.nan, np.nan, label=r'$\theta$', color='k', linewidth=2.2, linestyle=LiStyle)
                ax.plot(np.nan, np.nan, label=r'$\theta_{dagger}$', color='k', linewidth=0.8, linestyle=LiStyle)
            ax.legend()
            ax.grid(True)

            ###Subplot 2 = nothing for the moment
            axes[0, 1].axis('off')


            ###Subplot 3 : taub/N
            ax = axes[1,0]
            ###Plot horizontal line for value of C
            ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
            # Add text label above the line on the left side
            ax.text(
                35, C ,  # x and y coordinates in data units
                'C = {}'.format(C),
                ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                fontsize=10, color='gray'
            )
            ax.plot(interpolated_df['Time']*365, interpolated_df['taub']/interpolated_df['N'], color=colstar, linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.28, 0.57)
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
            ax = axes[2,0]
            ax.plot(interpolated_df['Time']*365, interpolated_df['ub']/365, color=colstar, linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.grid(True)


            ###Subplot 6 : Sigmann
            ax = axes[2,1]
            ax.plot(interpolated_df['Time']*365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.grid(True)

            # ax.plot(interpolated_df['Time']*365, interpolated_df['Chi'], color=colstar, linestyle=LiStyle, linewidth=2)
            # ax.axhspan(0, 2.2, color='grey', alpha=0.05)
            # ax.set_ylabel(r'$\chi$', fontsize=20)
            # ax.set_xlabel(r'Time [days]', fontsize=20)
            # ax.set_ylim(-3.5,2.2)
            # ax.grid(True)


            #### FIGURE 2 ####
            # Fill up figure 2 : Plot the ratio ub/N vs taub/N
            ax2.plot(interpolated_df['ub'] / (365*interpolated_df['N']), interpolated_df['taub'] / interpolated_df['N'],color=colstar,linestyle=LiStyle,linewidth=2)
            # Axis labels with LaTeX-style formatting
            ax2.set_xlabel(r'$\mathrm{u}_\mathrm{b} / N$ [m d$^{-1}$ MPa$^{-1}$]', fontsize=20)
            ax2.set_ylabel(r'$\tau_b / N$', fontsize=20)
            ax2.grid(True)

    ###Now add legend regarding C on empty upper right subplot
    axes[1, 1].legend(handles=handles,loc='center',fontsize=22,frameon=False)

    plt.show()