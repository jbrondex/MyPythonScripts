"""
@author: jbrondex

Description:
------------
This script generates a figure (for the supplementary) for the first version of the Nature paper led by Mylène Jacquemart.

The figure will likely consist of several subplots, showing map of various variable at bed over time

The exact layout and content of the figure will be refined as the analysis progresses.
"""
from cProfile import label

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
# Create legend handles manually, for example:
from matplotlib.lines import Line2D
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch
from matplotlib.colors import ListedColormap, BoundaryNorm

from pathlib import Path

import numpy as np
import re
import os
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
import matplotlib.dates as mdates
from datetime import datetime ###To handle dates

import matplotlib.ticker
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
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
#####################################################################################
# Define function to interpolate field defined on mesh nodes to the line profile ####
#####################################################################################
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y)
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

###############################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    #### SCRIPT PARAMETERS ###
    ### Parameter of Simulation of interest ###
    Cvalue= 0.3
    dcvalue = 5.0
    Year = '2025'
    StartSimuDay = 16 ####Day of May 2025 corresponding to simulation start
    SimuDuration = 15 ###Simu duration in days
    #EndOfLoadingDays = [19,20,21,22,23]
    EndOfLoadingDays = [23]
    ##in case we want to plot results for overload stopping before 23rd of May:
    ## One color per end of overloading day. In this order: blue/orange/green/brown/purple/grey/red/yellow
    ColorBlindFriendly_List = ['k'] ###['#377eb8', '#ff7f00', '#4daf4a', '#a65628', 'k'] ##'#984ea3','#999999', '#e41a1c', '#dede00']
    LineSizes = [2.4]###[1.8, 1.8, 1.8, 1.8, 2.4]

    ### Parameter of plot ###
    DaysOfMay25_ForMaps = [18, 20, 22, 24, 26] ##Day of May 2025 at which we want to plot the map of a variable
    ##Which maps do we want to show ? Can be 'ub', 'Sigma_nn', 'theta', 'taub_over_N','dtaub_over_dt','thetadag_minus_theta', 'taub_over_taud'
    variablesToMap = ['dtaub_over_dt', 'ub', 'theta', 'thetadag_minus_theta', 'taub_over_CN']###['ub', 'taub', 'Sigma_nn','theta', 'taub_over_N']##['ub', 'Sigma_nn']##Which maps do we want to show ? Can be 'ub', 'Sigma_nn', 'theta', 'taub_over_N'
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 30 ##in m/d !
    ###Do we want to show one/several star(s) on maps corresponding to position(s) of TimeLapse reconstructed velocities ?
    Stars = True ##True  ## Do we want to plot stars ?
    Stars_Size = 15 ###Size of star to show on maps
    Points_Timelapse_Ids = [3] ###Ids of points of TimeLapse reconstructed velocities that we want to plot (list with values among 1, 2, 3 and/or 4)
    valid_ids = {1, 2, 3, 4}
    Colors_Points_Timelapse = ['#cc00cc']##['#00bcd4']##['blue', 'orange','g','r'] ###Colors of dots and stars corresponding to the point at which TimeLapse velocity is reconstructed
    ### Contours we want to plot for debris thickness on maps
    DebrisContour_Levels = [10, 20, 30, 40, 50, 60]
    DebrisContour_LabelPositions = [(2630.73,1139.097), (2630.587,1139.089), (2630.877,1139.018), (2630.917,1138.978), (2630.736,1138.978), (2630.853,1138.930)]
    DebrisContour_Color = 'k'#'dimgray'#'lightcoral'
    DebrisContour_Thick = 0.9
    ###Defines bounds of the maps
    xmin,xmax = 2630250, 2631250
    ymin,ymax = 1138750, 1139500

    ########################
    ###TREAT OBSERVATIONS###
    ########################

    ####### TimeLapse Velocities ####
    ###Coordinates of surface points at which horizontal velocity is inferred
    Points_Timelapse = [[2630784, 1139171], [2630758, 1139194], [2630712, 1139233], [2630681, 1139269]]


    ###Load glacier contour
    Pathroot_Contour = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/1_Mesh/MESH{}/.'.format(Year))
    Filename_GlacierContour = 'Contour_Birch_{}_Julien.dat'.format(Year)
    Df_GlacierContour = pd.read_csv(Pathroot_Contour.joinpath(Filename_GlacierContour), names=['X', 'Y'], sep=r'\s+')
    xc = Df_GlacierContour['X'].values
    xc=np.append(xc,xc[0]) ##step required to close the contour
    yc = Df_GlacierContour['Y'].values
    yc = np.append(yc, yc[0])  ##step required to close the contour

    ### Load Orthophoto for Background
    # --- Load orthophoto
    Pathroot_Ortho = Path('/home/brondexj/Documents/Expertises/Blatten/OrthoPhotos/.')
    Filename_Ortho = 'OrthoPostCollapse.tif'
    with rasterio.open(Pathroot_Ortho.joinpath(Filename_Ortho)) as src:
        window = from_bounds(xmin, ymin, xmax, ymax, src.transform)
        ortho_crop = src.read(1, window=window)  # for grayscale
        # ortho_crop = src.read([1, 2, 3], window=window) # for RGB
        transform = src.window_transform(window)
    # --- Compute extent in km
    extent_km = [xmin / 1000, xmax / 1000, ymin / 1000, ymax / 1000]

    ### Load DEM of surface of glacier and rock deposit to build contour line of rock debris thickness
    Pathroot_DEMs= Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/1_Mesh/MESH2025/.')
    Filename_SurfGlacierDEM = 'DEM_surf_2025_Reconstructed_FromAdrien_MJacquemart_Smoothed.dat'
    Df_SurfGlacierDEM = pd.read_csv(Pathroot_DEMs.joinpath(Filename_SurfGlacierDEM), names=['X', 'Y', 'Z'], sep=r'\s+')
    Filename_SurfDebrisDEM = 'DEMSurfPre.dat'
    Df_SurfDebrisDEM = pd.read_csv(Pathroot_DEMs.joinpath(Filename_SurfDebrisDEM), names=['X', 'Y', 'Z'], sep=r'\s+')

    #### Load file containing factors to apply to Hdebris as a function of day of May for progressive load
    Pathroot_Factors= Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/.')
    Filename_Factors= 'OverloadFactors.dat'
    Df_Factors = pd.read_csv(Pathroot_Factors.joinpath(Filename_Factors), names=['Day','Factor'], sep=r'\s+', skiprows=1)

    #######################
    ###PREPARE THE PLOT ###
    ncols = len(DaysOfMay25_ForMaps)
    nrows = len(variablesToMap)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(3 * ncols, 3.5 * nrows),
        sharex='col',
        sharey=True,
        constrained_layout=True,
        squeeze=False
    )

    ########################
    ### START THE SCRIPT ###
    ####Rearrangement of parameters to get corresponding file names
    s = f"{Cvalue:.2f}"
    if len(s.rstrip('0').split('.')[1]) == 1:
        str_Cvalue = f"{int(Cvalue * 10):02d}"
    elif len(s.rstrip('0').split('.')[1]) == 2:
        str_Cvalue = f"{int(Cvalue * 100):03d}"

    if dcvalue <1:
        str_dcvalue = f"dc{int(dcvalue * 10):02d}m_"
    elif dcvalue.is_integer():
        str_dcvalue = f"dc{int(dcvalue):01d}m_"
    else:
        x, y =str(dcvalue).split(".")
        str_dcvalue = f"dc{x}dot{y}m_"

    ####definition of time limits from simu parameters
    t0 = pd.Timestamp("2025-05-{} 00:00".format(StartSimuDay))
    tf = t0 + pd.Timedelta(days=SimuDuration)
    CollapseTime = pd.Timestamp("2025-05-28 15:30")

    ####NOW LOAD DATA CORRESPONDING TO PRONOSTIC RUNS WITH OVERLOAD
    if Year == '2025':
        Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
        Filename = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue)
    else:
        Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat_{}/GlacierOut/.'.format(Year))
        Filename = 'MyBirch{}_RS_C{}_{}tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(Year, str_Cvalue,str_dcvalue)
    Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'Normal Stress', 'theta_dagger', 'theta', 'taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    ###Time step size for prono simu with overload:
    Timestep_size = 0.01 / 365  ### In years !!
    ###load file corresponding to first simu as dataframe
    print('Opening file:', Filename)
    Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
    ###Create a column with time (in years !)
    Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
    ###Convert absolute time in terms of date
    Df_wholebed['Date'] = t0 + pd.to_timedelta(Df_wholebed['Time'] * 365, unit='D')
    ###Divide dataframe in different df for the three boundaries
    Df_wholebed_front = Df_wholebed.loc[Df_wholebed['BC'] == 2].copy()
    Df_wholebed_bed = Df_wholebed.loc[Df_wholebed['BC'] == 3].copy()
    Df_wholebed_top = Df_wholebed.loc[Df_wholebed['BC'] == 4].copy()

    #####################################################################################################
    ###     PLOT MAPS OF VARIABLE OF INTEREST AT DAYS OF INTEREST FOR FIRST TWO LINES AND COLUMNS     ###
    #####################################################################################################
    #### For first and second lines we only work on dates of interest in May 25 (one day per column j)
    for j, DayOfMay25 in enumerate(DaysOfMay25_ForMaps ):
        ###Convert day in date
        DateOfInterest = pd.to_datetime("2025-05-{}".format(DayOfMay25))
        print('Dealing with maps for date:', DateOfInterest, 'May 2025')
        DateOfInterest_Prev = pd.to_datetime("2025-05-{}".format(DayOfMay25 - 1))
        ### Keep only lines for the date of interest
        ### For any given day of may of interest, I want to identify the hour that is present in the df and that is closest to the start of the day
        ExactDateOfInterest = Df_wholebed_bed[Df_wholebed_bed["Date"].dt.floor("D") == DateOfInterest]['Date'].min()
        ExactDateOfInterest_Prev = Df_wholebed_bed[Df_wholebed_bed["Date"].dt.floor("D") == DateOfInterest_Prev]['Date'].min()
        ### Keep only these lines
        Df_wholebed_SpecDay=Df_wholebed_bed[Df_wholebed_bed['Date']==ExactDateOfInterest].copy()
        Df_wholebed_SpecDay_Prev=Df_wholebed_bed[Df_wholebed_bed['Date']==ExactDateOfInterest_Prev].copy()
        ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes:
        Df_wholebed_SpecDay.loc[Df_wholebed_SpecDay['GM'] < 0, 'taub'] = 0
        Df_wholebed_SpecDay_Prev.loc[Df_wholebed_SpecDay_Prev['GM'] < 0, 'taub'] = 0
        ####RECOMPUTE BASAL VELOCITY####
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed_SpecDay['un'] = Df_wholebed_SpecDay['u'] * Df_wholebed_SpecDay['nx'] + Df_wholebed_SpecDay['v'] * Df_wholebed_SpecDay['ny'] + Df_wholebed_SpecDay['w'] * Df_wholebed_SpecDay['nz']
        ###Deduce ub from u and un
        Df_wholebed_SpecDay['ub'] = np.sqrt((Df_wholebed_SpecDay['u'] - Df_wholebed_SpecDay['un'] * Df_wholebed_SpecDay['nx']) ** 2 + (Df_wholebed_SpecDay['v'] - Df_wholebed_SpecDay['un'] * Df_wholebed_SpecDay['ny']) ** 2 + (Df_wholebed_SpecDay['w'] - Df_wholebed_SpecDay['un'] * Df_wholebed_SpecDay['nz']) ** 2)
        ####RECOMPUTE NORMAL STRESS FROM STOKES RESIDUAL ####
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed_SpecDay['Sigma_nn'] = (Df_wholebed_SpecDay['fx'] * Df_wholebed_SpecDay['nx'] + Df_wholebed_SpecDay['fy'] * Df_wholebed_SpecDay['ny'] +Df_wholebed_SpecDay['fz'] * Df_wholebed_SpecDay['nz']) / Df_wholebed_SpecDay['nodearea']
        ### Create an additionnal column that will be used to assess whether cavitation is at steady state or not (threshold at 98% pf steady)
        Df_wholebed_SpecDay['IsSteady'] = (abs(Df_wholebed_SpecDay['theta_dagger']-Df_wholebed_SpecDay['theta'])<0.02)
        ### Create an additionnal column to identify if we are in the coulomb regime locally (i.e. theta is steady AND taub = C*N with threshold of 0.02MPa = 20KPa)
        Df_wholebed_SpecDay['IsCoulomb'] = (Df_wholebed_SpecDay['IsSteady'] & (abs(Df_wholebed_SpecDay['taub']- Cvalue * (-Df_wholebed_SpecDay['Sigma_nn']))<0.02) & Df_wholebed_SpecDay['GM']>0)
        ###Create refined regular grid and interpolate field over grid (only once)
        if j ==0:
            xbed = np.arange(np.floor(np.min(Df_wholebed_SpecDay['X']))-10,np.floor(np.max(Df_wholebed_SpecDay['X']))+10,0.5)
            ybed = np.arange(np.floor(np.min(Df_wholebed_SpecDay['Y']))-10,np.floor(np.max(Df_wholebed_SpecDay['Y']))+10,0.5)
            Xbed, Ybed = np.meshgrid(xbed, ybed)
            ### Interpolate DEMs over same grid (just needed once)
            ZSurfGlacier = Interpolate_field(Df_SurfGlacierDEM, 'Z', Xbed, Ybed)
            ZSurfDebris =  Interpolate_field(Df_SurfDebrisDEM, 'Z', Xbed, Ybed)
            Hdebris_notclipped = ZSurfDebris - ZSurfGlacier
        ### For each considered day, get factor applied on debris thickness for progressive loading
        if DayOfMay25 < Df_Factors['Day'].min():
            f = 0
        elif DayOfMay25 > Df_Factors['Day'].max():
            f = Df_Factors[Df_Factors['Day'] == Df_Factors['Day'].max()]['Factor'].values
        else:
            f = Df_Factors[Df_Factors['Day'] == DayOfMay25]['Factor'].values
        Hdebris_fraction = f*Hdebris_notclipped
        ### If H debris is less than 1.5m -> Zero
        Hdebris = np.where(Hdebris_fraction < 1.5, 0, Hdebris_fraction)
        Hdebris_smoothed = gaussian_filter(Hdebris, sigma=6)
        ###Interpolate GM variable on refined regular grid as it is needed independently of variable to map
        GM = Interpolate_field(Df_wholebed_SpecDay,'GM',Xbed,Ybed)
        GM_mask = GM < 0
        GM_negative = (GM < 0).astype(float)
        ### From there we consider each variable to map independently (one variable per line i)
        ###~~~~~~~Loop on variables to plot~~~~~~~~~~~~~~~
        for i, VarToPlot in enumerate(variablesToMap):
            CS1 = None
            CS1_hatched = None
            CS1_GLcontour = None
            ###Get the proper subplot
            ax = axes[i, j]
            ###do some settings
            ax.tick_params(labelsize=9)  # fontsize of the tick labels
            ax.grid(True)
            ax.grid(alpha=0.5)
            ###Force x and y axis to same scale
            ax.set_aspect('equal', adjustable='box')
            ###Add title and label wherever needed
            if i == 0: ### First Line
                ax.set_title('{} May'.format(DayOfMay25), fontsize=21, weight='bold')
            elif i == len(variablesToMap): ###last line
                ax.set_xlabel(r'X [km]', fontsize=18)
            if j == 0:
                ax.set_ylabel(r'Y [km]', fontsize=18)
            # --- Plot orthophoto as background
            ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
            #####ADAPT THE PLOT DEPENDING ON VARIABLE THAT NEEDS TP BE PLOT
            if VarToPlot == 'ub':
                ###interpolate variable of interest on regular grid
                ub = Interpolate_field(Df_wholebed_SpecDay, 'ub', Xbed, Ybed)
                clevs=np.arange(0.0, velocitymax+0.1, 2.5)   ## cbar for shading
                Colormap = 'viridis'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res=5
                cbar_name = r'$\mathrm{u_b}$ [m d$^{-1}$]'
                CS1 = ax.contourf(Xbed/1000, Ybed/1000, ub/365, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
                CS1_GLcontour=None
            elif VarToPlot == 'taub':
                ###interpolate variable of interest on regular grid
                taub = Interpolate_field(Df_wholebed_SpecDay, 'taub', Xbed, Ybed)
                clevs=np.arange(0.0, 300.01, 50)   ## cbar for shading
                Colormap = 'magma_r'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res=50
                cbar_name = r'$\tau_\mathrm{b}$ [kPa]'
                VarToMask = taub*1000
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed/1000, Ybed/1000, VarMasked, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True,zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM_negative,levels=[0.5, 1.5],colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black',linewidths=0.8, linestyles='-',zorder=21)
                ### Highlight areas where taub is at upper bound
                #ax.scatter(Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsSteady']]['X'] / 1000,Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsSteady']]['Y'] / 1000,s=3,color='black',marker='.',zorder=10)
            elif VarToPlot == 'Sigma_nn':
                ###interpolate variable of interest on regular grid
                Sigma_nn = Interpolate_field(Df_wholebed_SpecDay, 'Sigma_nn', Xbed, Ybed)
                clevs = np.arange(-2.0, 2.01, 0.05)  ## cbar for shading
                vmin, vmax = clevs[0], clevs[-1]
                Colormap = 'RdBu'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res=0.5
                cbar_name = r'$\sigma_\mathrm{nn}$ [MPa]'
                VarToMask = Sigma_nn
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, VarMasked, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax, alpha=0.6,zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM_negative,levels=[0.5, 1.5],colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black',linewidths=0.8, linestyles='-',zorder=21)
                ### Highlight areas where taub is at upper bound
                #ax.scatter(Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['X'] / 1000,Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['Y'] / 1000,s=3,color='black',marker='.',zorder=10)
            elif VarToPlot == 'theta':
                ###interpolate variable of interest on regular grid
                theta = Interpolate_field(Df_wholebed_SpecDay,'theta',Xbed,Ybed)
                clevs = np.arange(0.0, 1.01, 0.1)  ## cbar for shading
                Colormap = 'cividis'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res=0.2
                cbar_name = r'$\theta$'
                VarToMask = 1 - theta
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, VarMasked, levels=clevs, cmap=cmap, extend='both', alpha=0.7,antialiased=True, zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM_negative,levels=[0.5, 1.5],colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black',linewidths=0.8, linestyles='-',zorder=21)
                ### Highlight areas where taub is at upper bound
                #ax.scatter(Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['X'] / 1000,Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['Y'] / 1000,s=3,color='black',marker='.',zorder=10)
            elif VarToPlot == 'taub_over_CN':
                ###interpolate variable of interest on regular grid
                taub = Interpolate_field(Df_wholebed_SpecDay, 'taub', Xbed, Ybed)
                Sigma_nn = Interpolate_field(Df_wholebed_SpecDay, 'Sigma_nn', Xbed, Ybed)
                ### Make the division by being careful of places where N might be zero or NaN
                taub_over_CN_raw = np.divide(taub, -Sigma_nn*Cvalue, out=np.full_like(taub, np.nan), where=(Sigma_nn != 0))
                # Mask NaNs and infs for safe plotting
                taub_over_CN = np.ma.masked_invalid(taub_over_CN_raw)
                clevs=np.arange(0.0, 2.01, 0.25)
                Colormap = 'RdBu'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res=0.5
                cbar_name = r'$\tau_\mathrm{b}/ \mathrm{C N}$'
                CS1 = ax.contourf(Xbed/1000, Ybed/1000, taub_over_CN, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
                # # White + hatched filling where GM < 0
                # CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM,levels=[-1.0, -0.1],colors='lightgrey', zorder=5)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black',linewidths=0.8, linestyles='-',zorder=21)
            elif VarToPlot == 'dtaub_over_dt':
                ###interpolate variable of interest on regular grid
                taub = Interpolate_field(Df_wholebed_SpecDay, 'taub', Xbed, Ybed)
                taub_Prev = Interpolate_field(Df_wholebed_SpecDay_Prev, 'taub', Xbed, Ybed)
                clevs=np.arange(-100, 101, 10.0)
                Colormap = 'coolwarm'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res=50.0
                cbar_name = r'$\Delta \tau_\mathrm{b}$ [kPa]'
                VarToMask = (taub-taub_Prev)*1000
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed/1000, Ybed/1000, VarToMask, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True, zorder=2)
                # # White + hatched filling where GM < 0
                # CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM_negative,levels=[0.5, 1.5],colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black',linewidths=0.8, linestyles='-',zorder=21)
            elif VarToPlot == 'thetadag_minus_theta':
                ###interpolate variable of interest on regular grid
                theta = Interpolate_field(Df_wholebed_SpecDay,'theta',Xbed,Ybed)
                theta_dagger = Interpolate_field(Df_wholebed_SpecDay,'theta_dagger',Xbed,Ybed)
                ###Here I make my own colormap to show clear transition from steady state areas (i.e. thetadag-theta<0.02) from transient areas
                clevs = [
                    0.0, 0.005, 0.01, 0.015, 0.02,
                    0.05, 0.1, 0.25, 0.5, 0.75, 1.0
                ]
                colors = [
                    '#08306b',
                    '#2171b5',
                    '#6baed6',
                    '#c6dbef',  # < 0.01 : blues

                    '#ffffcc',
                    '#ffeda0',
                    '#fed976',
                    '#feb24c',
                    '#fd8d3c',
                    '#e31a1c'  # > 0.01 : yellow -> red
                ]
                cmap = ListedColormap(colors)
                norm = BoundaryNorm(clevs, cmap.N)
                cbar_ticks=[0, 0.01, 0.02, 0.1, 0.5, 1.0]
                cbar_name = r'$\theta^{*} - \theta$'
                VarToMask = ((1-theta_dagger) -(1 - theta))
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, VarMasked , levels=clevs, cmap=cmap, norm = norm, extend='both', alpha=0.7,antialiased=True, zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM_negative,levels=[0.5, 1.5],colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black',linewidths=0.8, linestyles='-',zorder=21)
            else:
                raise ValueError(f"Invalid value in variablesToMap: '{VarToPlot}'. "
                        "Expected one of: 'ub', 'Sigma_nn', 'theta', 'taub_over_CN','dtaub_over_dt','thetadag_minus_theta', 'taub_over_taud'.")
            ###NOW DO WHAT IS INDEPENDENT OF VARIABLE TO PLOT ###
            vmin, vmax = clevs[0], clevs[-1]
            ### Plot contours of debris thickness
            CS1_debriscontour = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed, levels=DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-', zorder=10)
            ###Put labels with position determined manually
            # ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d", manual=DebrisContour_LabelPositions, fontsize=7.5)
            ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d", fontsize=6.5)
            ###Glacier contour
            ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2, zorder=22)
            ###Plot stars where we have reconstructed velocity from TimeLapse camera on ub map
            if Stars and VarToPlot == 'ub':
                for Points_Timelapse_Id, colstar in zip(Points_Timelapse_Ids, Colors_Points_Timelapse):
                    if len(Points_Timelapse_Ids) != len(Colors_Points_Timelapse):
                        raise ValueError(
                            "Points_Timelapse_Ids and Colors_Points_Timelapse "
                            "must have the same length."
                        )
                    if Points_Timelapse_Id not in valid_ids:
                        raise ValueError(f"Invalid value in Points_Timelapse_Ids: '{Points_Timelapse_Id}'. "
                        "There are 4 points at which we have data so should be 1, 2, 3 or 4.")
                    xstar, ystar = Points_Timelapse[Points_Timelapse_Id-1]
                    ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Stars_Size,markeredgecolor='black', markeredgewidth=1.2)
            ###Below we remove colors that are outside of glacier contour
            clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
            patch = PathPatch(clippath, facecolor='none')
            ax.add_patch(patch)
            for c in [CS1, CS1_GLcontour, CS1_hatched, CS1_debriscontour]:
                if c is not None:
                    c.set_clip_path(patch)
            if j == len(DaysOfMay25_ForMaps) - 1:
                cbar = fig.colorbar(CS1, ax=ax, shrink=0.98, pad=0.05)
                if VarToPlot == 'thetadag_minus_theta':
                    cbar.set_ticks(cbar_ticks)
                else:
                    cbar.set_ticks(np.arange(vmin, vmax+0.001, cbar_ticks_Res))
                cbar.ax.tick_params(labelsize=16)
                cbar.set_label(cbar_name, fontsize=21)
    plt.show()