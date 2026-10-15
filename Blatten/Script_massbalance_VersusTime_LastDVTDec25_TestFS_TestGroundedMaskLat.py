"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mass balance over time for chosen C values with and without FreeSurface solved
and with or without the groundedmask lateral activated
It compares volume variation to flux at boundaires and FS residuals when FS are solved

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
from matplotlib.lines import Line2D

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
    # Cases = ['NoFS NoGMLat','NoFS GMLat', 'FS NoGMLat', 'FS GMLat']
    # LineStyles = ['-.',':','--','-']
    Cases = ['NoFS NoGMLat','NoFS GMLat', 'FS NoGMLat', 'FS GMLat']
    LineStyles = ['-.',':','--','-']
    LineSize = 2.4 ###Default value of linewidth
    ###When FS is solved what value of alpha for stabilization for the upper FS
    # Alphas = ['10m1','10m2']
    # LineSizes = [2.4, 1.5]
    Alphas = ['10m2']
    LineSizes = [2]
    ###PREPARE THE PLOT
    legend_handles = []
    ###~~~~~~~FIG 1 and 1 bis (one per C value) : six subplots of mean flux across the four boundaries + residuals at top and bed [m3/m2/jour]~~~~~~~~~~~~~~~
    ###Prepare the subplot
    ###Prepare the figure : Mean flux at boundaries : lateral, front, top, bottom
    fig1, axes1 = plt.subplots(nrows=3, ncols=2, figsize=(15, 12), sharex=True)
    axes_flat1 = axes1.flatten()
    fig1bis, axes1bis = plt.subplots(nrows=3, ncols=2, figsize=(15, 12), sharex=True)
    axes_flat1bis = axes1bis.flatten()
    ###~~~~~~~FIG 2 : Relative volume over time versus total flux + residuals ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig2, axes2 = plt.subplots(nrows=2, ncols=2, figsize=(15, 12), sharex=True)  # , gridspec_kw={'hspace': -0.15})
    axes_flat2 = axes2.flatten()


    ### START OF SCRIPT###
    ###Loop on each C considered
    for i, (C, blue_col) in enumerate(zip(Cvalues, blues)):
        print('Dealing with C value :', C)
        str_Cvalue = f"{int(C * 10):02d}"
        ###Loop on Case considered (NoFS versus FS)
        for j, (Case, LiStyle) in enumerate(zip(Cases, LineStyles)):
            print('Dealing with case :', Case)
            for k, (Alpha, LineSize) in enumerate(zip(Alphas, LineSizes)):
                if Case in ['NoFS NoGMLat', 'NoFS GMLat']: ###No alpha if no FS solved
                    if k > 0:
                        continue
                else:
                    print('Dealing with Alpha :', Alpha)
                if Case == 'NoFS NoGMLat':
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_LASTDVT_Dec2025/2_RateAndState_NoFS_NoGMLat/GlacierOut/.')
                    Filename_mb = 'MyBirch_RS_C{}_NoFS_NoGMLat_tsp001d__masscons.dat'.format(str_Cvalue)
                    ###For the save scalar file used for the mass balance closure
                    Col_Name_mb = ['Timestep', 'Volume', 'FluxLat', 'FluxFront', 'FluxBed', 'FluxSurf', 'MinFlux','MaxFlux', 'AreaLat', 'AreaFront', 'AreaBed', 'AreaSurf']
                    # pattern = re.compile(r"MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(str_Cvalue))
                elif Case=='NoFS GMLat' and k==0:
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_LASTDVT_Dec2025/3_RateAndState_NoFS_GMLat/GlacierOut/.')
                    Filename_mb = 'MyBirch_RS_C{}_NoFS_GMLat_tsp001d__masscons.dat'.format(str_Cvalue)
                    ###For the save scalar file used for the mass balance closure
                    Col_Name_mb = ['Timestep', 'Volume', 'FluxLat', 'FluxFront', 'FluxBed', 'FluxSurf', 'MinFlux','MaxFlux', 'AreaLat', 'AreaFront', 'AreaBed', 'AreaSurf']
                elif Case=='FS NoGMLat':
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_LASTDVT_Dec2025/4_RateAndState_FS_NoGMLat/GlacierOut/.')
                    ###For the save scalar file used for the mass balance closure
                    Col_Name_mb = ['Timestep', 'Volume', 'FluxLat', 'FluxFront', 'FluxBed', 'FluxSurf', 'MinFlux','MaxFlux', 'AreaLat', 'AreaFront', 'AreaBed', 'AreaSurf', 'ZsResidual', 'ZbResidual']
                    Filename_mb = 'MyBirch_RS_C{}_BothFS_NoGMLat_tsp001d_StabFabFSTop_Alpha{}_FSBotStabStdy__masscons.dat'.format(str_Cvalue,Alpha)
                elif Case=='FS GMLat':
                    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_LASTDVT_Dec2025/5_RateAndState_FS_GMLat/GlacierOut/.')
                    ###For the save scalar file used for the mass balance closure
                    Col_Name_mb = ['Timestep', 'Volume', 'FluxLat', 'FluxFront', 'FluxBed', 'FluxSurf', 'MinFlux','MaxFlux', 'AreaLat', 'AreaFront', 'AreaBed', 'AreaSurf', 'ZsResidual', 'ZbResidual']
                    Filename_mb = 'MyBirch_RS_C{}_BothFS_GMLat_tsp001d_StabFabFSTop_Alpha{}_FSBotStabStdy__masscons.dat'.format(str_Cvalue,Alpha)
                ###Correct time step size for the simu for which I changes to 0.01 day:
                Timestep_size = 0.01 / 365  ### In years !!
                ###load file as dataframe
                Df_wholebed_mb = pd.read_csv(Pathroot.joinpath(Filename_mb), names=Col_Name_mb, sep=r'\s+')
                ###########################################################################
                ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED C AND CONSIDERED CASE####
                ###Create a column with time (in years !)
                Df_wholebed_mb['Time'] = Timestep_size * Df_wholebed_mb['Timestep']
                ###/!\/!\ BE VERY CAREFUL: FLUX ENTERING THE DOMAIN ARE DEFINED AS NEGATIVE IN ELMER (because negative normal)
                ###/!\/!\ HERE I SWITCH THIS CONVENTION TO HAVE A SOUND MASS BALANCE: POSITIVE FLUX = MASS GAIN, NEGATIVE FLUX = MASS LOSS
                BCFlux=['FluxLat', 'FluxFront', 'FluxBed', 'FluxSurf']
                Df_wholebed_mb[BCFlux] *= -1
                ###Calculate mean flux across each of the boundaries
                for BC in ['Lat', 'Front', 'Bed', 'Surf']:
                    Df_wholebed_mb[f'Flux{BC}_perArea'] = (Df_wholebed_mb[f'Flux{BC}'] / Df_wholebed_mb[f'Area{BC}'])
                    ###also calculate the cumulated flux at each boundary (i.e. integrated overtime)
                    Df_wholebed_mb[f'CumFlux{BC}'] = (Df_wholebed_mb[f'Flux{BC}'] * Timestep_size).cumsum()
                ###For the with FS case, calculate mean residuals across top and bottom boundaries
                if 'ZsResidual' in Col_Name_mb:
                    Df_wholebed_mb['ZsResidual_perArea'] = (Df_wholebed_mb['ZsResidual'] / Df_wholebed_mb['AreaSurf'])
                    ###also calculate the cumulated residual  (i.e. integrated overtime)
                    Df_wholebed_mb[f'CumZsResidual'] = (Df_wholebed_mb['ZsResidual'] * Timestep_size).cumsum()
                if 'ZbResidual' in Col_Name_mb:
                    Df_wholebed_mb['ZbResidual_perArea'] = (Df_wholebed_mb['ZbResidual'] / Df_wholebed_mb['AreaBed'])
                    ###also calculate the cumulated residual  (i.e. integrated overtime)
                    Df_wholebed_mb[f'CumZbResidual'] = (Df_wholebed_mb['ZbResidual'] * Timestep_size).cumsum()
                ##############################################################
                ######### NOW WHAT DO WE WANT TO PLOT ????  ##################
                ##############################################################
                #### Two plots (one per C) containing 6 subplots with flux at BC (Lat, Front, Bed, Surf, Bed Residuals, Surf Residuals)
                #### For all considered cases
                for k1,Flux in enumerate(['FluxLat', 'FluxFront', 'FluxBed', 'FluxSurf', 'ZbResidual', 'ZsResidual']):
                    ###Residual are meaningfull only for case with FS
                    if 'Residual' in Flux and 'NoFS' in Case:
                        continue
                    if C == 0.3:
                        ###Get proper ax
                        ax = axes_flat1[k1]
                        ### After 5 days of simu with C =0.3 the simu are in the choux
                        ax.set_xlim(0, 5.01)
                        # Major ticks at specific values
                        ax.set_xticks([1, 2, 3, 4, 5])
                        # Minor ticks every 1 unit
                        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
                        ###ylim depends on C value and Flux considered
                        if Flux == 'FluxLat' or Flux == 'FluxFront':
                            ax.set_ylim(-6, 6)
                        elif Flux == 'FluxBed' or Flux == 'FluxSurf':
                            ax.set_ylim(-2, 2)
                        elif Flux == 'ZbResidual' or Flux == 'ZsResidual':
                            ax.set_ylim(-0.5, 0.5)
                    elif C == 0.5:
                        ###Get proper ax
                        ax = axes_flat1bis[k1]
                        ### For C =0.5 the simu are generally running until the end
                        ax.set_xlim(0, 24.01)
                        # Major ticks at specific values
                        ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
                        # Minor ticks every 1 unit
                        ax.xaxis.set_minor_locator(MultipleLocator(1))
                        ###ylim depends on C value and Flux considered
                        if Flux == 'FluxLat' or Flux == 'FluxFront':
                            ax.set_ylim(-4, 4)
                        elif Flux == 'FluxBed' or Flux == 'FluxSurf':
                            ax.set_ylim(-1, 1)
                        elif Flux == 'ZbResidual' or Flux == 'ZsResidual':
                            ax.set_ylim(-0.1, 0.1)
                    #### Other settings are common to all C
                    ax.tick_params(labelsize=14)  # fontsize of the tick labels
                    ax.grid(True)
                    ax.grid(alpha=0.5)
                    ax.set_title('{}'.format(Flux), fontsize=21, weight='bold')
                    # Add grid for both major and minor ticks
                    ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                    ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)
                    if k1 in [0, 2, 4]:
                        ax.set_ylabel(r'Mean Flux [$\mathrm{m/d}$]', fontsize=20)
                    if k1 in [4,5]:
                        ax.set_xlabel(r'Time [days]', fontsize=20)
                    ###Plot the flux in question versus time
                    ax.plot(Df_wholebed_mb['Time']*365, Df_wholebed_mb[f'{Flux}_perArea']/365, color=blue_col, linestyle=LiStyle, linewidth=LineSize)
                    #### Prepare a common legend
                    if k1 == 0 and i == 0:
                        legend_handles.append(Line2D([0], [0],color='k',linestyle=LiStyle,linewidth=LineSize,label=Case))

                ################################################################################################
                ### FIG2 : PLOT showing volume over time for all combi tested (one per subplot both C on all subplots)
                ################################################################################################
                ###Get proper ax
                ax = axes_flat2[j]
                #### Other settings are common to all C
                ax.tick_params(labelsize=14)  # fontsize of the tick labels
                ax.grid(True)
                ax.grid(alpha=0.5)
                ax.set_title('{}'.format(Case), fontsize=21, weight='bold')
                ax.set_xlim(0, 24.01)
                ax.set_ylim(-0.5, 0.5)
                # Major ticks at specific values
                ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
                # Minor ticks every 1 unit
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)
                if j in [0, 2]:
                    ax.set_ylabel(r'$\Delta$V [$\mathrm{Mm^3}$]', fontsize=20)
                if j in [2,3]:
                    ax.set_xlabel(r'Time [days]', fontsize=20)
                ###Plot the flux in question versus time
                ax.plot(Df_wholebed_mb['Time'] * 365, (Df_wholebed_mb['Volume'] - Df_wholebed_mb[Df_wholebed_mb['Timestep']==1]['Volume'].values) / 10**6,color=blue_col,linestyle='-', linewidth=LineSize)
                ax.plot(Df_wholebed_mb['Time'] * 365, Df_wholebed_mb['CumFluxFront'] / 10**6, color=blue_col,linestyle='-', linewidth=0.8)
                if Case=='FS NoGMLat' or Case=='FS GMLat':
                    ax.plot(Df_wholebed_mb['Time'] * 365, (Df_wholebed_mb['CumFluxLat'] + Df_wholebed_mb['CumFluxFront'] + Df_wholebed_mb['CumFluxBed'] + Df_wholebed_mb['CumFluxSurf'] + Df_wholebed_mb['CumZsResidual'] + Df_wholebed_mb['CumZbResidual']) / 10 ** 6, color=blue_col, linestyle='--',linewidth=LineSize)
                else:
                    ax.plot(Df_wholebed_mb['Time'] * 365, (Df_wholebed_mb['CumFluxLat'] + Df_wholebed_mb['CumFluxFront'] + Df_wholebed_mb['CumFluxBed'] + Df_wholebed_mb['CumFluxSurf']  ) / 10 ** 6, color=blue_col, linestyle='--',linewidth=LineSize)
                ###make one common legend
                if i ==0:
                    legend_handles2 = [
                    Line2D([0], [0], color='k', linestyle='-',  linewidth=LineSize, label='Relative Volume'),
                    Line2D([0], [0], color='k', linestyle='-',  linewidth=0.8,      label='Cum. Flux at front'),
                    Line2D([0], [0], color='k', linestyle='--', linewidth=LineSize, label=r'$\Sigma$ Cum. Fluxes'),
                    ]

#### --- BUILD FINAL LEGEND WITH CORRECT ORDER ---
    fig1.legend(handles=legend_handles,fontsize=18,loc='upper center', ncol=len(legend_handles), bbox_to_anchor=(0.5, 1.0))
    fig1bis.legend(handles=legend_handles,fontsize=18,loc='upper center', ncol=len(legend_handles), bbox_to_anchor=(0.5, 1.0))
    fig2.legend(handles=legend_handles2,fontsize=18,loc='upper center', ncol=len(legend_handles), bbox_to_anchor=(0.5, 1.0))
    #             ###Calculate ub from velocity and normal vector
    #             ###Start by calculating un (should be zero but there is some residual)
    #             Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
    #             ###Deduce ub from u and un
    #             Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
    #             ### Calculate mean ub over whole glacier for each time step
    #             mean_ub_by_time = Df_wholebed.groupby('Time')['ub'].mean()
    #             ### Show mean ub in Fig 1
    #             if LiStyle == '-' and k==0:
    #                 ax.plot(mean_ub_by_time.index*365, mean_ub_by_time/365, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=LineSize)
    #             else:
    #                 ax.plot(mean_ub_by_time.index * 365, mean_ub_by_time / 365, color=blue_col, linestyle=LiStyle, linewidth=LineSize)
    #             ###Dummy plots for legend regarding lateral BC
    #             if i==0:
    #                 if Case=='FS GMLat':
    #                     ax.plot(np.nan,np.nan,label='Case = {}, Alpha ={}'.format(Case, Alpha),color='k',linestyle=LiStyle, linewidth=LineSize)
    #                 else:
    #                     ax.plot(np.nan,np.nan,label='Case = {}'.format(Case),color='k',linestyle=LiStyle, linewidth=LineSize)
    #             ax.set_ylim(0.0, 40.0)
    #             ax.set_xlim(0, 24.01)
    #             # Major ticks at specific values
    #             ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
    #             # Minor ticks every 1 unit
    #             ax.xaxis.set_minor_locator(MultipleLocator(1))
    #             # Add grid for both major and minor ticks
    #             ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
    #             ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)
    #             ###Plot evolution of grounded area over time in Fig2
    #             ### Calculate total area over the whole glacier for each time step
    #             Total_Area = Df_wholebed.groupby('Time')['nodearea'].sum()
    #             ### Grounded area
    #             Total_GroundedArea = (Df_wholebed[Df_wholebed['GM'] == 1].groupby('Time')['nodearea'].sum())
    #             ### Fraction of grounded area over total area
    #             Ratio = Total_GroundedArea/Total_Area
    #             ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
    #             # Copy ratio but mask values where Ratio == 1
    #             Ratio_masked = Ratio.copy()
    #             Ratio_masked[Ratio_masked == 1] = np.nan
    #
    #             ### Plot ratio
    #             if LiStyle == '-':
    #                 ax2.plot(Total_Area.index*365, Ratio_masked, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=2.4)
    #             else:
    #                 ax2.plot(Total_Area.index*365, Ratio_masked, color=blue_col,linestyle=LiStyle, linewidth=2.4)
    #             ###Dummy plots for legend regarding Case
    #             if i==0:
    #                 ax2.plot(np.nan,np.nan,label='Case = {}'.format(Case),color='k',linestyle=LiStyle, linewidth=2.4)
    #             ax2.set_xlim(0, 24.01)
    #             # Major ticks at specific values
    #             ax2.set_xticks([0, 4, 8, 12, 16, 20, 24])
    #             # Minor ticks every 1 unit
    #             ax2.xaxis.set_minor_locator(MultipleLocator(1))
    #             # Add grid for both major and minor ticks
    #             ax2.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
    #             ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)
    # # --- BUILD FINAL LEGEND WITH CORRECT ORDER ---
    # handles, labels = ax.get_legend_handles_labels()
    # # Separate C legends (colored solid lines) from BC legends (black dummy)
    # handles_C = []
    # labels_C = []
    # handles_BC = []
    # labels_BC = []
    # handles_other = []
    # labels_other = []
    #
    # for h, l in zip(handles, labels):
    #     if l.startswith("C ="):
    #         handles_C.append(h)
    #         labels_C.append(l)
    #     elif l.startswith("Case ="):
    #         handles_BC.append(h)
    #         labels_BC.append(l)
    #     else:
    #         handles_other.append(h)
    #         labels_other.append(l)
    #
    # # Combine them: C first, BC after
    # handles_ordered = handles_C + handles_BC + handles_other
    # labels_ordered = labels_C + labels_BC + labels_other
    #
    # ax.legend(handles_ordered, labels_ordered, fontsize=22)
    # ax2.legend(handles_ordered, labels_ordered, fontsize=22, loc='lower right')
    #


    plt.show()