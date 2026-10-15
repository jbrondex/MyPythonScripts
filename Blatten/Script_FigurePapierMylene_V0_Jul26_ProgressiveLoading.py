"""
@author: jbrondex

Description:
------------
This script generates the figures for the first version of the Nature paper led by Mylène Jacquemart.

The figure will likely consist of several subplots, including:
- two maps of ice velocity at different dates,
- two maps of either τ_b/N or σ_nn at different dates,
- the temporal evolution of the modeled front velocity compared with radar observations,
- possibly the temporal evolution of modeled velocities at Ethan's measurement points compared with
  velocities reconstructed from the timelapse imagery.

The exact layout and content of the figure will be refined as the analysis progresses.
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
# Create legend handles manually, for example:
from matplotlib.lines import Line2D

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

    ### Parameter of plot ###
    ##Parameters for maps of first subplots
    DaysOfMay25_ForMaps = [21, 24] ##Day of May 2025 at which we want to plot the map of a variable
    ##Which maps do we want to show ? Can be 'ub', 'Sigma_nn', 'theta', 'taub_over_CN','dtaub_over_dt','thetadag_minus_theta', 'taub_over_taud'
    variablesToMap = ['dtaub_over_dt','ub']##
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 30 ##in m/d !
    ###Do we want to show one/several star(s) on maps corresponding to position(s) of TimeLapse reconstructed velocities ?
    Stars = True ##True  ## Do we want to plot stars ?
    Stars_Size = 15 ###Size of star to show on maps
    Points_Timelapse_Ids = [3] ###Ids of points of TimeLapse reconstructed velocities that we want to plot (list with values among 1, 2, 3 and/or 4)
    valid_ids = {1, 2, 3, 4}
    Colors_Points_Timelapse = ['#cc00cc']##['#00bcd4']##['blue', 'orange','g','r'] ###Colors of dots and stars corresponding to the point at which TimeLapse velocity is reconstructed
    ###Do we want to show hatches in places where taub = CN and theta is at steady state ?
    IsHatches = False
    ### Contours we want to plot for debris thickness on maps
    DebrisContour_Levels = [10, 20, 30, 40, 50, 60]
    DebrisContour_LabelPositions = [(2630.73,1139.097), (2630.587,1139.089), (2630.877,1139.018), (2630.917,1138.978), (2630.736,1138.978), (2630.853,1138.930)]
    DebrisContour_Color = 'k'#'dimgray'#'lightcoral'
    DebrisContour_Thick = 0.9
    ###Defines bounds of the maps
    xmin,xmax = 2630250, 2631250
    ymin,ymax = 1138750, 1139500
    ###For the last subplot, do we want velocities in log scale ?
    LogScale = False #True

    ##Parameters last subplot (velocity versus time)
    ###For bottom plot, chose color and linesize of full overload case
    Color_FullOverloadCase = 'k'
    LineSize_FullOverloadCase = 2.4
    ###Limit on yaxis for last subplot
    ylimlastsubplot = 60
    #### Do we want to plot the cases for which overload is stopped on a given day ?
    PlotPartialOverload = False #True
    EndOfLoadingDays = [19, 20, 21,22]  ###Which cases of partial overload do we want to plot ? Full overload is systematically plotted
    ## One color per end of overloading day. In this order: blue/orange/green/brown/purple/grey/red/yellow
    ColorBlindFriendly_List = ['#377eb8', '#ff7f00', '#4daf4a', '#a65628']
    LineSizes = [1.8, 1.8, 1.8, 1.8]

    ########################
    ###TREAT OBSERVATIONS###
    ######## Radar Data #########
    RadarPosition = [2628689, 1141335,1794] ##xradar, yradar, zradar
    LastReliableRadarData = pd.Timestamp("2025-05-25 06:00") ###Radar data is not reliable anymore after 6 am on the 25th of May 25
    Col_Name_Radar = ['Time','Fullfront_6h','Middlefront_6h','Upperfront_6h','Fullfront_24h','Middlefront_24h','Upperfront_24h','Fullfront_RWB','Middlefront_RWB','Upperfront_RWB','FullMean']
    Pathroot_Radar = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/Obs/.')
    Filename_Radar = 'glacier_front_RadarData_mmperday.csv'
    ###BE WARE : Radar data are in mm/day
    Df_Radar = pd.read_csv(Pathroot_Radar.joinpath(Filename_Radar), header=0, names=Col_Name_Radar)
    ###Convert time of csv file in panda date
    Df_Radar['Date'] = pd.to_datetime(Df_Radar["Time"], format="%d.%m.%y %H:%M")
    ###Cut all data that are not reliable (i.e. aquired later than last reliable data)
    Df_Radar =  Df_Radar[Df_Radar['Date']<LastReliableRadarData]
    ###Calculate a mean between the three processing method (6h vs 24h vs RWB) for the full front
    Df_Radar["mean_Fullfront"] = Df_Radar[["Fullfront_6h", "Fullfront_24h", "Fullfront_RWB"]].mean(axis=1, skipna=True)

    ####### TimeLapse Velocities ####
    ###Coordinates of surface points at which horizontal velocity is inferred
    Points_Timelapse = [[2630784, 1139171], [2630758, 1139194], [2630712, 1139233], [2630681, 1139269]]
    Col_Name_Timelapse = ['X0','Y0','t0','t1','horizontal_velo','sigma']
    Pathroot_Timelapse = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/Obs/.')
    Filename_Timelapse = 'VelocitiesFromTimelapse.csv'
    Df_Timelapse = pd.read_csv(Pathroot_Timelapse.joinpath(Filename_Timelapse), header=0, names=Col_Name_Timelapse)
    ###Get the mean time between two images used to compute velocities
    Df_Timelapse['t0'] = pd.to_datetime(Df_Timelapse['t0'], utc=True)
    Df_Timelapse['t1'] = pd.to_datetime(Df_Timelapse['t1'], utc=True)
    Df_Timelapse['t_mid'] = Df_Timelapse['t0'] + (Df_Timelapse['t1'] - Df_Timelapse['t0']) / 2

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
    #### Convert May day numbers into actual dates for plot of stepwise overload increase
    DebrisDates = pd.to_datetime([f'2025-05-{day:02d}' for day in Df_Factors['Day']])
    ###Add the last day of simu in DebrisDates
    EndSimuDate = pd.Timestamp(year=2025, month=5, day=StartSimuDay)+ pd.Timedelta(days=SimuDuration)
    DebrisDates = DebrisDates.append(pd.DatetimeIndex([EndSimuDate]))
    ###Now deal with overload factors
    OverloadFactors=Df_Factors['Factor'].values
    OverloadFactors = np.append(OverloadFactors, OverloadFactors[-1])
    TotalDebrisVolume=3.36 ###Total debris volume on 23rd of May in Mm3
    #######################
    ###PREPARE THE PLOT ###
    fig = plt.figure(figsize=(10, 12), constrained_layout=True)
    gs = fig.add_gridspec(nrows=3,ncols=2,height_ratios=[1, 1, 1.2])
    # Top row
    ax00 = fig.add_subplot(gs[0, 0])
    ax01 = fig.add_subplot(gs[0, 1], sharey=ax00)
    # Second row
    ax10 = fig.add_subplot(gs[1, 0], sharex=ax00, sharey=ax00)
    ax11 = fig.add_subplot(gs[1, 1], sharex=ax01, sharey=ax10)
    # Store map axes in a 2D list
    axes = [
        [ax00, ax01],
        [ax10, ax11]
    ]
    # # Bottom row spanning both columns
    ax20 = fig.add_subplot(gs[2, :])
    # Right column: hide y tick labels
    for ax in (ax01, ax11):
        plt.setp(ax.get_yticklabels(), visible=False)
    # Top row: hide x tick labels
    for ax in (ax00, ax01):
        plt.setp(ax.get_xticklabels(), visible=False)

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
    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
    Filename = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad23May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue)
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
        Df_wholebed_SpecDay['IsSteady'] = (abs(Df_wholebed_SpecDay['theta_dagger']-Df_wholebed_SpecDay['theta'])<0.03)
        ### If wanted, create an additionnal column to identify if we are in the coulomb regime locally (i.e. theta is steady AND taub = C*N with threshold of 0.02MPa = 20KPa)
        if IsHatches:
            ratio = (Df_wholebed_SpecDay['taub']/ (Cvalue * -Df_wholebed_SpecDay['Sigma_nn']))
            Df_wholebed_SpecDay['IsCoulomb'] = (Df_wholebed_SpecDay['IsSteady'] & (ratio.between(0.92, 1.08)) & (Df_wholebed_SpecDay['GM']>0))
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
        ### Interpolate IsCoulomb over same grid (needed just once)
        if IsHatches:
            IsCoulomb = griddata((Df_wholebed_SpecDay['X'], Df_wholebed_SpecDay['Y']),Df_wholebed_SpecDay['IsCoulomb'].astype(float), (Xbed, Ybed), method='nearest')
        ### From there we consider each variable to map independently (one variable per line i)
        ###~~~~~~~Loop on variables to plot~~~~~~~~~~~~~~~
        for i, VarToPlot in enumerate(variablesToMap):
            CS1 = None
            CS1_Coulomb = None
            CS1_hatched = None
            CS1_GLcontour = None
            CS1_debriscontour = None
            ###Get the proper subplot
            ax = axes[i][j]
            ###do some settings
            ax.tick_params(labelsize=12)  # fontsize of the tick labels
            ax.grid(True)
            ax.grid(alpha=0.5)
            ###Force x and y axis to same scale
            ax.set_aspect('equal', adjustable='box')
            ###Add title and label wherever needed
            if i == 0: ### First Line
                ax.set_title('{} May'.format(DayOfMay25), fontsize=21, weight='bold')
            else: ###second line
                ax.set_xlabel(r'X [km]', fontsize=18)
            if j == 0:
                ax.set_ylabel(r'Y [km]', fontsize=18)
            # --- Plot orthophoto as background
            ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
            #####ADAPT THE PLOT DEPENDING ON VARIABLE THAT NEEDS TP BE PLOT
            if VarToPlot == 'ub':
                ###interpolate variable of interest on regular grid
                ub = Interpolate_field(Df_wholebed_SpecDay, 'ub', Xbed, Ybed)
                clevs = np.arange(0.0, velocitymax + 0.1, 2.5)  ## cbar for shading
                Colormap = 'viridis'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res = 5
                cbar_name = r'$\mathrm{u_b}$ [m d$^{-1}$]'
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, ub / 365, levels=clevs, cmap=cmap, extend='both', alpha=0.7,antialiased=True)
                CS1_GLcontour = None
            elif VarToPlot == 'taub':
                ###interpolate variable of interest on regular grid
                taub = Interpolate_field(Df_wholebed_SpecDay, 'taub', Xbed, Ybed)
                clevs = np.arange(0.0, 300.01, 50)  ## cbar for shading
                Colormap = 'magma_r'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res = 50
                cbar_name = r'$\tau_\mathrm{b}$ [kPa]'
                VarToMask = taub * 1000
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, VarMasked, levels=clevs, cmap=cmap, extend='both',alpha=0.7, antialiased=True, zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000, Ybed / 1000, GM_negative, levels=[0.5, 1.5], colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black', linewidths=0.8,linestyles='-', zorder=21)
                ### Highlight areas where taub is at upper bound
                # ax.scatter(Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsSteady']]['X'] / 1000,Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsSteady']]['Y'] / 1000,s=3,color='black',marker='.',zorder=10)
            elif VarToPlot == 'Sigma_nn':
                ###interpolate variable of interest on regular grid
                Sigma_nn = Interpolate_field(Df_wholebed_SpecDay, 'Sigma_nn', Xbed, Ybed)
                clevs = np.arange(-2.0, 2.01, 0.05)  ## cbar for shading
                vmin, vmax = clevs[0], clevs[-1]
                Colormap = 'RdBu'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res = 0.5
                cbar_name = r'$\sigma_\mathrm{nn}$ [MPa]'
                VarToMask = Sigma_nn
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, VarMasked, cmap=cmap, shading='auto', vmin=vmin,vmax=vmax, alpha=0.6, zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000, Ybed / 1000, GM_negative, levels=[0.5, 1.5], colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black', linewidths=0.8,linestyles='-', zorder=21)
                ### Highlight areas where taub is at upper bound
                # ax.scatter(Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['X'] / 1000,Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['Y'] / 1000,s=3,color='black',marker='.',zorder=10)
            elif VarToPlot == 'theta':
                ###interpolate variable of interest on regular grid
                theta = Interpolate_field(Df_wholebed_SpecDay, 'theta', Xbed, Ybed)
                clevs = np.arange(0.0, 1.01, 0.1)  ## cbar for shading
                Colormap = 'cividis'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res = 0.2
                cbar_name = r'$\theta$'
                VarToMask = 1 - theta
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, VarMasked, levels=clevs, cmap=cmap, extend='both',alpha=0.7, antialiased=True, zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000, Ybed / 1000, GM_negative, levels=[0.5, 1.5], colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black', linewidths=0.8,linestyles='-', zorder=21)
                ### Highlight areas where taub is at upper bound
                # ax.scatter(Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['X'] / 1000,Df_wholebed_SpecDay[Df_wholebed_SpecDay['IsCoulomb']]['Y'] / 1000,s=3,color='black',marker='.',zorder=10)
            elif VarToPlot == 'taub_over_CN':
                ###interpolate variable of interest on regular grid
                taub = Interpolate_field(Df_wholebed_SpecDay, 'taub', Xbed, Ybed)
                Sigma_nn = Interpolate_field(Df_wholebed_SpecDay, 'Sigma_nn', Xbed, Ybed)
                ### Make the division by being careful of places where N might be zero or NaN
                taub_over_CN_raw = np.divide(taub, -Sigma_nn * Cvalue, out=np.full_like(taub, np.nan),where=(Sigma_nn != 0))
                # Mask NaNs and infs for safe plotting
                taub_over_CN = np.ma.masked_invalid(taub_over_CN_raw)
                clevs = np.arange(0.0, 2.01, 0.25)
                Colormap = 'RdBu'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res = 0.5
                cbar_name = r'$\tau_\mathrm{b}/ \mathrm{C N}$'
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, taub_over_CN, levels=clevs, cmap=cmap, extend='both',alpha=0.7, antialiased=True)
                # # White + hatched filling where GM < 0
                # CS1_hatched = ax.contourf(Xbed / 1000,Ybed / 1000,GM,levels=[-1.0, -0.1],colors='lightgrey', zorder=5)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black', linewidths=0.8,linestyles='-', zorder=21)
                ## Hatch areas where basal friction is in Coulomb regime
                if IsHatches:
                    plt.rcParams['hatch.linewidth'] = 1.4###thicker hatches
                    CS1_Coulomb = ax.contourf(Xbed / 1000,Ybed / 1000,IsCoulomb,levels=[0.5, 1.5],colors='none', hatches=['////'], zorder=20)
            elif VarToPlot == 'dtaub_over_dt':
                ###interpolate variable of interest on regular grid
                taub = Interpolate_field(Df_wholebed_SpecDay, 'taub', Xbed, Ybed)
                taub_Prev = Interpolate_field(Df_wholebed_SpecDay_Prev, 'taub', Xbed, Ybed)
                clevs = np.arange(-100, 101, 10.0)
                Colormap = 'coolwarm'
                cmap = cmx.get_cmap(Colormap)
                cbar_ticks_Res = 50.0
                cbar_name = r'$\Delta \tau_\mathrm{b}$ [kPa]'
                VarToMask = (taub - taub_Prev) * 1000
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, VarToMask, levels=clevs, cmap=cmap, extend='both',alpha=0.7, antialiased=True, zorder=2)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black', linewidths=0.8,linestyles='-', zorder=21)
                ## Hatch areas where basal friction is in Coulomb regime
                if IsHatches:
                    plt.rcParams['hatch.linewidth'] = 1.4###thicker hatches
                    CS1_Coulomb = ax.contourf(Xbed / 1000,Ybed / 1000,IsCoulomb,levels=[0.5, 1.5],colors='none', hatches=['////'], zorder=20)
            elif VarToPlot == 'thetadag_minus_theta':
                ###interpolate variable of interest on regular grid
                theta = Interpolate_field(Df_wholebed_SpecDay, 'theta', Xbed, Ybed)
                theta_dagger = Interpolate_field(Df_wholebed_SpecDay, 'theta_dagger', Xbed, Ybed)
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
                cbar_ticks = [0, 0.01, 0.02, 0.1, 0.5, 1.0]
                cbar_name = r'$\theta^{*} - \theta$'
                VarToMask = ((1 - theta_dagger) - (1 - theta))
                VarMasked = np.ma.masked_where(GM_mask, VarToMask)
                CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, VarMasked, levels=clevs, cmap=cmap, norm=norm,extend='both', alpha=0.7, antialiased=True, zorder=2)
                # White + hatched filling where GM < 0
                CS1_hatched = ax.contourf(Xbed / 1000, Ybed / 1000, GM_negative, levels=[0.5, 1.5], colors='lightgrey',zorder=20)
                ### At floating areas (GM <0), make a contour
                CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.01], colors='black', linewidths=0.8,linestyles='-', zorder=21)
            else:
                raise ValueError(f"Invalid value in variablesToMap: '{VarToPlot}'. "
                                 "Expected one of: 'ub', 'Sigma_nn', 'theta', 'taub_over_CN','dtaub_over_dt','thetadag_minus_theta', 'taub_over_taud'.")
            ###NOW DO WHAT IS INDEPENDENT OF VARIABLE TO PLOT ###
            vmin, vmax = clevs[0], clevs[-1]
            if VarToPlot == 'ub':
                ### Plot contours of debris thickness
                CS1_debriscontour = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed, levels=DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-')
                ###Put labels with position determined manually
                ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d", fontsize=6)
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
            for c in [CS1, CS1_GLcontour, CS1_Coulomb, CS1_hatched, CS1_debriscontour]:
                if c is not None:
                    c.set_clip_path(patch)
            if j == len(DaysOfMay25_ForMaps) - 1:
                cbar = fig.colorbar(CS1, ax=ax, shrink=0.95, pad=0.05)
                if VarToPlot == 'thetadag_minus_theta':
                    cbar.set_ticks(cbar_ticks)
                else:
                    cbar.set_ticks(np.arange(vmin, vmax+0.001, cbar_ticks_Res))
                cbar.ax.tick_params(labelsize=14)
                cbar.set_label(cbar_name, fontsize=20)

    #########################################################################
    ###     PLOT EVOLUTION OF VELOCITY VS TIME IN TWO BOTTOM SUBPLOTS     ###
    #########################################################################
    ####FRONT VELOCITY IN RADAR DIRECTION ####
    ###for each node of front, we calculate the vector from radar to node:
    dx = Df_wholebed_front['X'] - RadarPosition[0]
    dy = Df_wholebed_front['Y'] - RadarPosition[1]
    dz = Df_wholebed_front['Z'] - RadarPosition[2]
    ###for each node of front, we calculate norm of this vector:
    r_norm = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    ###Now do the projection through scalar product with unit vector
    Df_wholebed_front['uf_radarproj'] = (Df_wholebed_front['u'] * dx + Df_wholebed_front['v'] * dy + Df_wholebed_front['w'] * dz) / r_norm
    ###Now calculate the mean of all nodes at front for each Date
    mean_uf_radarproj_by_time = Df_wholebed_front.groupby('Date')['uf_radarproj'].mean()
    ####SURFACE VELOCITY####
    ###first calculate mean surface velocity
    Df_wholebed_top['us_norm'] = np.sqrt(Df_wholebed_top['u'] ** 2 + Df_wholebed_top['v'] ** 2 + Df_wholebed_top['w'] ** 2)
    ###Now calculate the mean norm of surface velocity for each Date
    mean_us_norm_by_time = Df_wholebed_top.groupby('Date')['us_norm'].mean()
    ###Now I interpolate the modeled horizontal velocity at the points of the timelapse camera for each date
    results = []
    #first I will do that for each date available
    for t, group in Df_wholebed_top.groupby('Date'):
        ###get x and y of all nodes
        nodes = group[['X', 'Y']].values
        # get horizontal velocities at nodes
        u_vals = group['u'].values
        v_vals = group['v'].values
        # do the interpolation at each desired oints of the timelapse camera
        for Points_Timelapse_Id in Points_Timelapse_Ids:
            if Points_Timelapse_Id not in valid_ids:
                raise ValueError(f"Invalid value in Points_Timelapse_Ids: '{Points_Timelapse_Id}'. "
                "There are 4 points at which we have data so should be 1, 2, 3 or 4.")
            Point = Points_Timelapse[Points_Timelapse_Id-1]
            u_interp = griddata(nodes, u_vals, Point, method='linear')[0]
            v_interp = griddata(nodes, v_vals, Point, method='linear')[0]
            results.append({
                'Date': t,
                'Point_id': Points_Timelapse_Id, ##id starts at 1
                'Xpt': Point[0],
                'Ypt': Point[1],
                'u_interp': u_interp,
                'v_interp': v_interp
            })
    Df_ModelInterpAtTimeLapsePoints = pd.DataFrame(results)
    ###Calculate norm of horizontal velocity
    Df_ModelInterpAtTimeLapsePoints['u_interp_norm'] = np.sqrt(Df_ModelInterpAtTimeLapsePoints['u_interp'] ** 2 + Df_ModelInterpAtTimeLapsePoints['v_interp'] ** 2)

    ### ~~~~~ From there we can start the subplot ~~~~~~~
    ###Get the proper subplot
    ax = ax20
    ax.tick_params(labelsize=14)
    ax.set_ylabel(r'Ice velocity [m d$^{-1}$]', fontsize=18)
    ### If desired, use a logarithmic scale for all velocities
    if LogScale:
        ax.set_yscale('log')
    ###We plot Mean front velo in radar direction + Horizontal velocity at selected surface points ~~~~~~~~~~~~~~~
    ### Radar Data
    ax.plot(Df_Radar['Date'], Df_Radar['mean_Fullfront'].abs()/1000, color='b', label="Radar observations (front)", linestyle='--', linewidth=LineSize_FullOverloadCase*1.5, zorder=20)
    ### Modeled mean front velo in radar direction
    ax.plot(mean_uf_radarproj_by_time.index, mean_uf_radarproj_by_time.abs()/365, label="Modelled front velocity", color='k',linestyle='--', linewidth=LineSize_FullOverloadCase)
    ### TimeLapse Velocities Modeled versus Observed
    for Points_Timelapse_Id, colstar in zip(Points_Timelapse_Ids, Colors_Points_Timelapse):
        if len(Points_Timelapse_Ids) != len(Colors_Points_Timelapse):
            raise ValueError(
                "Points_Timelapse_Ids and Colors_Points_Timelapse "
                "must have the same length."
            )
        if Points_Timelapse_Id not in valid_ids:
            raise ValueError(f"Invalid value in Points_Timelapse_Ids: '{Points_Timelapse_Id}'. "
                             "There are 4 points at which we have data so should be 1, 2, 3 or 4.")
        Point = Points_Timelapse[Points_Timelapse_Id-1]
        subset2 = Df_ModelInterpAtTimeLapsePoints[Df_ModelInterpAtTimeLapsePoints['Point_id'] == Points_Timelapse_Id]
        subset_Obs = Df_Timelapse[(Df_Timelapse['X0'] == Point[0]) & (Df_Timelapse['Y0'] == Point[1])]
        ax.scatter(subset_Obs['t_mid'],subset_Obs['horizontal_velo'], label="Time-lapse observations (surface point)", color=colstar, s=20)
        ax.plot(subset2['Date'],subset2['u_interp_norm'] / 365, label="Modelled velocity (surface point)", color='k',linestyle='-',linewidth=LineSize_FullOverloadCase)
    ###Add vertical line for collapse day and annotate
    ax.axvline(CollapseTime, linestyle='--', linewidth=2.7, color='r')
    ax.text(CollapseTime+ pd.Timedelta(hours=17, minutes=30), ylimlastsubplot-5,'Collapse', fontsize=18, color='r', rotation=90, va='top', ha='right')
    ###Add vertical lines for the dates shown on the two columns of subplots above
    # ax.axvline(pd.to_datetime("2025-05-{}".format(DaysOfMay25_ForMaps[0])), color='0.5', ls=':', lw=1, zorder=10)
    # ax.axvline(pd.to_datetime("2025-05-{}".format(DaysOfMay25_ForMaps[1])), color='0.5', ls=':', lw=1, zorder=10)
    # #ax.text(CollapseTime, 40,'Effondrement', fontsize=26, color='r', rotation=90, va='top', ha='right')
    ###Dummy plots for legend regarding lateral B
    if LogScale:
        ax.set_ylim(0.01, 300.0)
    else:
        ax.set_ylim(0.0, ylimlastsubplot)
    ax.set_xlim(t0, tf)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
    if LogScale:
        ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10))
        ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10,subs=np.arange(2, 10) * 0.1))
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    # Add grid for both major and minor ticks
    ax.grid(which="major", color="black", linestyle="-", linewidth=0.9, alpha=0.6)
    ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7, alpha=0.4)
    for label in ax20.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')

    ####Here I deal with the right axis about debris overload
    # Secondary y-axis
    ax_debris = ax.twinx()
    # Step-like debris loading
    ax_debris.fill_between(DebrisDates,OverloadFactors*TotalDebrisVolume, step='post', alpha=0.18, color='grey', zorder=0, linewidth=1.0)
    ax_debris.step(DebrisDates,OverloadFactors*TotalDebrisVolume,where='post',color='grey',linewidth=1.0,zorder=1)
    ax_debris.set_ylabel('Cumulative debris\nvolume [Mm$^3$]',fontsize=18,color='0.5')
    ax_debris.tick_params(axis='y',colors='0.5',labelsize=14)
    ax_debris.set_ylim(0, 4.0)

    ax.set_zorder(ax_debris.get_zorder() + 1)
    ax.patch.set_visible(False)

    ####NOW I WILL PLOT THE PARTIAL OVERLOAD CASES IF REQUIRED
    if PlotPartialOverload:
        ###Loop on each C considered
        for k, (EndDay, blue_col, LineSize) in enumerate(zip(EndOfLoadingDays, ColorBlindFriendly_List, LineSizes)):
            ####NOW LOAD DATA CORRESPONDING TO PRONOSTIC RUNS WITH OVERLOAD
            print('Dealing with simu with overload finished on :', EndDay, 'of May')
            Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Juin2026/4_PRONO_PROGRESSIVELOAD_FS_NoGMLat/GlacierOut/.')
            Filename_PartialOv = 'MyBirch_RS_C{}_{}tsp001d_Start16May_StopLoad{}May_Stop31May__bed.dat'.format(str_Cvalue,str_dcvalue,EndDay)
            ###load file corresponding to first simu as dataframe
            print('Opening file:', Filename_PartialOv )
            Df_wholebed_PartialOv = pd.read_csv(Pathroot.joinpath(Filename_PartialOv), names=Col_Name, sep=r'\s+')
            ###########################################################################
            ###NOW WE CAN WORK ON THE DATAFRAME OF CONSIDERED END OF OVERLOAD DAY
            ###Create a column with time (in years !)
            Df_wholebed_PartialOv['Time'] = Timestep_size * Df_wholebed_PartialOv['Timestep']
            ###Convert absolute time in terms of date
            Df_wholebed_PartialOv['Date'] = t0 + pd.to_timedelta(Df_wholebed_PartialOv['Time'] * 365, unit='D')
            ###Create a dedicated dataframe for each boundary of interest
            Df_wholebed_PartialOv_front = Df_wholebed_PartialOv[Df_wholebed_PartialOv['BC'] == 2]
            Df_wholebed_PartialOv_bed = Df_wholebed_PartialOv[Df_wholebed_PartialOv['BC'] == 3]
            Df_wholebed_PartialOv_top = Df_wholebed_PartialOv[Df_wholebed_PartialOv['BC'] == 4]

            ####SURFACE VELOCITY####
            ###first calculate mean surface velo
            Df_wholebed_PartialOv_top['us_norm'] = np.sqrt(Df_wholebed_PartialOv_top['u'] ** 2 + Df_wholebed_PartialOv_top['v'] ** 2 + Df_wholebed_PartialOv_top['w'] ** 2)
            ###Now calculate the mean norm of surface velocity for each Date
            mean_us_norm_by_time_PartialOv = Df_wholebed_PartialOv_top.groupby('Date')['us_norm'].mean()
            ###Now I interpolate the modeled horizontal velocity at the four points of the timelapse camera for each date
            results = []
            # first I will do that for each date available
            for t, group in Df_wholebed_PartialOv_top.groupby('Date'):
                ###get x and y of all nodes
                nodes = group[['X', 'Y']].values
                # get horizontal velocities at nodes
                u_vals = group['u'].values
                v_vals = group['v'].values
                # do the interpolation at each of the four points of the timelapse camera
                for id, Point in enumerate(Points_Timelapse):
                    u_interp = griddata(nodes, u_vals, Point, method='linear')[0]
                    v_interp = griddata(nodes, v_vals, Point, method='linear')[0]
                    results.append({
                        'Date': t,
                        'Point_id': id + 1,  ##id starts at 1
                        'Xpt': Point[0],
                        'Ypt': Point[1],
                        'u_interp': u_interp,
                        'v_interp': v_interp
                    })
            Df_ModelInterpAtTimeLapsePoints_PartialOv = pd.DataFrame(results)
            ###Calculate norm of horizontal velocity
            Df_ModelInterpAtTimeLapsePoints_PartialOv['u_interp_norm'] = np.sqrt(Df_ModelInterpAtTimeLapsePoints_PartialOv['u_interp'] ** 2 + Df_ModelInterpAtTimeLapsePoints_PartialOv['v_interp'] ** 2)
            ### TimeLapse Velocities Modeled versus Observed
            for Points_Timelapse_Id, colstar in zip(Points_Timelapse_Ids, Colors_Points_Timelapse):
                if len(Points_Timelapse_Ids) != len(Colors_Points_Timelapse):
                    raise ValueError(
                        "Points_Timelapse_Ids and Colors_Points_Timelapse "
                        "must have the same length."
                    )
                if Points_Timelapse_Id not in valid_ids:
                    raise ValueError(f"Invalid value in Points_Timelapse_Ids: '{Points_Timelapse_Id}'. "
                                     "There are 4 points at which we have data so should be 1, 2, 3 or 4.")
                Point = Points_Timelapse[Points_Timelapse_Id - 1]
                subset2_PartialOv = Df_ModelInterpAtTimeLapsePoints_PartialOv[Df_ModelInterpAtTimeLapsePoints_PartialOv['Point_id'] == Points_Timelapse_Id]
                ax.plot(subset2_PartialOv['Date'],subset2_PartialOv['u_interp_norm'] / 365, label="Modelled velocity (surface point)", color=blue_col,linestyle='-',linewidth=LineSize)

    ###############################################
    ### NOW DEAL WITH LEGEND FOR BOTTOM SUBPLOT ###
    ###############################################
    if PlotPartialOverload:
        legend_handles = [
            Line2D([0], [0],color='k',linestyle='-',linewidth=LineSize_FullOverloadCase,label='Modelled velocity (surface point)'),
            Line2D([0], [0],color=colstar,marker='o',linestyle='None',markersize=5,label='Observed velocity (surface point)'),
            Line2D([0], [0],color='b',linestyle='--',linewidth=LineSize_FullOverloadCase*1.5,label='Observed front velocity (radar)'),
            Line2D([0], [0],color='k',linestyle='--',linewidth=LineSize_FullOverloadCase,label='Modelled mean front velocity'),
            # Empty object used as a section title
            Line2D([], [],linestyle='None',label='Overload stopped on:'),
            Line2D([0], [0],color=ColorBlindFriendly_List[0],linestyle='-',linewidth=LineSizes[0],label='19 May'),
            Line2D([0], [0],color=ColorBlindFriendly_List[1],linestyle='-',linewidth=LineSizes[1],label='20 May'),
            Line2D([0], [0],color=ColorBlindFriendly_List[2],linestyle='-',linewidth=LineSizes[2],label='21 May'),
            Line2D([0], [0],color=ColorBlindFriendly_List[3],linestyle='-',linewidth=LineSizes[3],label='22 May')
        ]
    else:
        legend_handles = [
            Line2D([0], [0],color='k',linestyle='-',linewidth=LineSize_FullOverloadCase,label='Modelled velocity (surface point)'),
            Line2D([0], [0],color=colstar,marker='o',linestyle='None',markersize=5,label='Observed velocity (surface point)'),
            Line2D([0], [0],color='b',linestyle='--',linewidth=LineSize_FullOverloadCase*1.5,label='Observed front velocity (radar)'),
            Line2D([0], [0],color='k',linestyle='--',linewidth=LineSize_FullOverloadCase,label='Modelled mean front velocity')
            ]
    leg = ax20.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.98),
        fontsize=9,
        frameon=True,
        borderaxespad=0
    )
    
    ###To have a line in bold
    for text in leg.get_texts():
        if text.get_text() == 'Overload stopped on:':
            text.set_fontweight('bold')



    plt.show()