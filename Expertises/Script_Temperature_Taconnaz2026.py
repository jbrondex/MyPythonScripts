"""
@author: jbrondex

Description:
------------
This file aims at : 1/analyse the transient evolution of temperatures over the 2026 summer campaign
2/ Show continuous temperature field at each borehole over the whole season
3/ Find steady state temperature and compare then to that measured during previous campaigns


"""
from json.decoder import NaN
from math import nan

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import matplotlib.dates as mdates


from pathlib import Path

import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
from datetime import datetime, date, timedelta ###To handle dates
import matplotlib.dates as mdates
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
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

###This function return the slope of T=ax+b where x is a rolling window of size 7 days by default
def rolling_slope(series, window=7):
    x = np.arange(window)

    return series.rolling(window).apply(
        lambda y: np.polyfit(x, y, 1)[0],
        raw=True
    )

if __name__ == "__main__":
    #### Parameters of script ####
    Forage2026_ID = ['1', '10', '3', '4']
    Forage2026_Depth = [80.8, 83.4, 103, 89.10]
    Forage2026_XYLambert2 = [[949532.598,2105271.283],[949432.983,2105118.072],[949231.325,2104511.397],[949558.043,2104002.718]]
    Forage2026_AltIGN69_RAF20=[3411.45,3491.22,3785.36,4064.10]
    Forage2026_RelPosThermFromBottomToTop = [[0,5,10,15,20,25,30,40,50,60,70,80,90,100,110],[0,5,10,15,20,25,30,40,50,60,70,80,90,100,110],[0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100],[0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100]]
    Tlim_ColorBar =[-2.5, -2.5, -7, -10]

    LimitDepth = -2 ###We do not consider any thermistors shallower than this limit
    XdayWindow = 4 ###Size in days of rolling window over which to calculate T decrease

    # Colors for past temperature profiles
    YearColors = {
        2008: 'blue',
        2012: 'green',
        2013: 'cyan',
        2014: 'orange',
        2021: 'purple',
        2026: 'red'
    }

    ###Do we want to plot a map of interpolated basal temperature
    InterpBasalT = True
    if InterpBasalT:
        ###Do we want to have 2 independent figures for 2021 and 2026 or a subplot with both ?
        PlotAsSubplots = True
        ###limit of map
        xmin = 948700
        xmax = 950800
        ymin = 2103500
        ymax= 2105700
        ### Load Orthophoto for Background
        # --- Load orthophoto
        Pathroot_Ortho = Path('/home/brondexj/Documents/Expertises/Taconnaz/CampagneTerrain2026/ProjetQgis_Taco2026/.')
        Filename_Ortho = 'TacoMap.tif'
        with rasterio.open(Pathroot_Ortho.joinpath(Filename_Ortho)) as src:
            window = from_bounds(xmin, ymin, xmax, ymax, src.transform)
            # ortho_crop = src.read(1, window=window)  # for grayscale
            ortho_crop = src.read([1, 2, 3], window=window) # for RGB
            ortho_crop = np.moveaxis(ortho_crop, 0, -1)
            transform = src.window_transform(window)
        # --- Compute extent in km
        extent_km = [xmin / 1000, xmax / 1000, ymin / 1000, ymax / 1000]
        ###Color associated to each borehole on map
        BoreholeColors = {
            'Forage_01': '#e31a1c',
            'Forage_03': '#33a02c',
            'Forage_04': '#2c0000',
            'Forage_06': '#38cee3',
            'Forage_07': '#800080',
            'Forage_08': '#ff8c00',
            'Forage_09': '#ffbb1c',
            'Forage_10': '#ff00ff',
            '25-H': '#954aa2',
        }
        BoreholeLabels = {
            'Forage_01': '01',
            'Forage_03': '03',
            'Forage_04': '04',
            'Forage_06': '06',
            'Forage_07': '07',
            'Forage_08': '08',
            'Forage_09': '09',
            'Forage_10': '10',
            '25-H': '25-H',
        }

    #### START OF SCRIPT####
    ###Set stuff for steady T that do not need to be re-evaluated for each borehole
    ColNames_2026_SteadyT = ['Borehole','Year','X','Y','Alt_Surf','Borehole_Depth','Alt_Bed','Sensor_ID', 'Alt_Therm','Sensor_Depth', 'T']
    ColNames_Past_SteadyT = ['Borehole','Year','X','Ywrong','Y','Alt_Surf','Sensor_Depth','T']
    Pathroot_Steady = Path('/home/brondexj/Documents/Expertises/Taconnaz/CampagneTerrain2026/Données2026/SteadyTemperatures/.')
    FileName_Past_SteadyT = 'Forages_2008-2021_all.csv'
    ###Open file for past campaign that gather dat for all boreholes
    Df_Past_SteadyT = pd.read_csv(Pathroot_Steady.joinpath(FileName_Past_SteadyT), skiprows=1, header=None, names=ColNames_Past_SteadyT,sep=',')

    ThermDepth_2026 = {}
    if InterpBasalT:
        DeepestSensors2026 = []
        Df_2021_SteadyT = Df_Past_SteadyT[Df_Past_SteadyT['Year'] == 2021]
        Df_DeepestSensors2021 = Df_2021_SteadyT.loc[Df_2021_SteadyT.groupby('Borehole')['Sensor_Depth'].idxmin()].reset_index(drop=True)
        ###Remove bis where there are any
        Df_DeepestSensors2021['Borehole'] = (Df_DeepestSensors2021['Borehole'].str.replace('bis', '', regex=False))
    for i,(ID,Depth,XY,Alt,RelPos, Tlim_ColBar) in enumerate(zip(Forage2026_ID,Forage2026_Depth,Forage2026_XYLambert2, Forage2026_AltIGN69_RAF20,Forage2026_RelPosThermFromBottomToTop, Tlim_ColorBar)):
        #####################################################
        ### START WITH MEASUREMENT FILES AT 15MN TIMESTEP ###
        #####################################################
        if ID in ['1']:
            ColNames_2026 = ['Time','Number','Batt_Volt','PTemp', 'T15', 'T14', 'T13', 'T12', 'T11', 'T10', 'T9', 'T8', 'T7', 'T6', 'T5', 'T4', 'T3', 'T2', 'T1']
        elif ID in ['10']:
            ColNames_2026 = ['Time','Number','Batt_Volt','PTemp', 'T15', 'T14', 'T13', 'T12', 'T11', 'T10', 'T9', 'T8', 'T7', 'T6', 'T5', 'T4', 'T3', 'T2', 'T1', 'Tair']
        elif ID in ['3','4']:
            ColNames_2026 = ['Time','Number','Batt_Volt','PTemp', 'T16', 'T15', 'T14', 'T13', 'T12', 'T11', 'T10', 'T9', 'T8', 'T7', 'T6', 'T5', 'T4', 'T3', 'T2', 'T1']
        # Depth of each thermistor
        ThermDepth_2026[ID] = {f'T{j + 1}': -Depth + RelPos[j] for j in range(len(RelPos))}
        # Keep only thermistors that are sufficiently deep
        ThermBelowSurface = [therm for therm, depth in ThermDepth_2026[ID].items() if depth < LimitDepth]
        Pathroot = Path('/home/brondexj/Documents/Expertises/Taconnaz/CampagneTerrain2026/Données2026/Thermistances_01092026/Forage_{}/.'.format(ID))
        FileNames_2026 = '20260903_Taco-T{}_temp_15min.dat'.format(ID)
        ###Open file
        Df_FullTherm_2026 = pd.read_csv(Pathroot.joinpath(FileNames_2026), skiprows=4, header=None, names=ColNames_2026, sep=',',parse_dates=['Time'], na_values=['NAN'])
        ### Keep only column of interest, i.e. time and thermistors that are sufficiently deep
        Df_Therm_2026 = Df_FullTherm_2026[['Time'] + ThermBelowSurface].copy()

        #######################################
        ### PLOT dT/dt to find steady state ###
        #######################################
        ### Calculate daily mean temperatures
        Df_Daily = (Df_Therm_2026.set_index('Time').resample('1D').mean())
        ### From there calculate the mean slope of T=f(t) (in °C/day) over a rolling window of X days centered on the considered day
        Df_SlopeXd = Df_Daily.apply(rolling_slope, window=XdayWindow)
        ###Plot the evolution of mean slope for each thermistors to find steady state
        fig, axes = plt.subplots(len(Df_SlopeXd.columns), 1,figsize=(12, 2.5 * len(Df_SlopeXd.columns)),sharex=True)
        for ax, therm in zip(axes, reversed(Df_SlopeXd.columns)):
            ax.axhline(0.0, linestyle='-', color='k', linewidth=0.75)
            # Thermistor depth
            depth = ThermDepth_2026[ID][therm]
            ax.plot(Df_SlopeXd.index,Df_SlopeXd[therm],label=f'{depth:.0f} m')
            # Stationarity thresholds
            ax.axhline(0.01, linestyle='--', color='k', linewidth=0.75)
            ax.axhline(-0.01, linestyle='--', color='k', linewidth=0.75)
            # ax.set_ylim([-0.2, 0.2])
            ax.legend()
            ax.grid(alpha=0.3)
        # Common labels
        fig.supylabel(r'Slope [$^\circ$C d$^{-1}$]', fontsize=12)
        axes[-1].set_xlabel('Date')
        axes[-1].xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO, interval=1))
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
        # Borehole title
        fig.suptitle(f'Borehole {ID}', fontsize=16)
        plt.setp(axes[-1].get_xticklabels(), rotation=45, ha='right')
        plt.tight_layout(rect=[0.01, 0.05, 1, 0.96])

        #######################################################################################
        ### PLOT T=f(t) continuously over full borehole depth and full measurement campaign ###
        #######################################################################################
        ####Replace outliers (thermistors broken) by NaN in Df_Therm_2026
        Df_Therm_clean = Df_Therm_2026.copy()
        # Resample daily average at a 15mn timestep for comparison to full timeserie
        Daily_mean_15min = Df_Daily.reindex(Df_Therm_2026['Time'], method='ffill')
        # Give it the same index as Df_Therm_2026
        Daily_mean_15min.index = Df_Therm_2026.index
        # Create a mask for outliers (T above 1°C, T below -12°C, more than 2°C gap between T daily and T15mn)
        for therm in ThermBelowSurface:
            mask_bad = ((Df_Therm_2026[therm] > 1) |(Df_Therm_2026[therm] < -12) |(abs(Df_Therm_2026[therm] - Daily_mean_15min[therm]) > 3))
            # Replace by NaN
            Df_Therm_clean.loc[mask_bad, therm] = np.nan
        # ---------------------------------------------------------
        # Get thermistors and corresponding depths
        # ---------------------------------------------------------
        # Thermistors below the surface
        therms = ThermBelowSurface
        # Corresponding depths
        depths = np.array([ThermDepth_2026[ID][therm]for therm in therms])
        # Sort from deepest to shallowest
        order = np.argsort(depths)
        depths = depths[order]
        therms = np.array(therms)[order]
        # ---------------------------------------------------------
        # Time and temperature arrays
        # ---------------------------------------------------------
        dates = Df_Therm_clean['Time']
        T = Df_Therm_clean[therms].to_numpy()
        # ---------------------------------------------------------
        # Regular vertical grid for interpolation
        # ---------------------------------------------------------
        depth_grid = np.linspace(-Depth,0,300)
        # Array that will contain interpolated temperatures
        # dimensions = depth x time
        T_interp = np.full((len(depth_grid), len(dates)),np.nan)
        # ---------------------------------------------------------
        # Interpolate temperature vertically at each timestep
        # ---------------------------------------------------------
        for i in range(len(dates)):
            temp = T[i, :]
            # Remove possible NaN values
            valid = np.isfinite(temp)
            if valid.sum() >= 2:
                f = interp1d(depths[valid],temp[valid],kind='linear',bounds_error=False,fill_value=np.nan)
                T_interp[:, i] = f(depth_grid)
        # ---------------------------------------------------------
        # Plot
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(14, 7))
        cmap = plt.get_cmap('rainbow').copy()
        # Couleur des valeurs masquées (NaN)
        cmap.set_bad(color='lightgrey')
        im = ax.pcolormesh(dates,depth_grid,T_interp,shading='auto',cmap=cmap, vmin=Tlim_ColBar, vmax=0)
        # ---------------------------------------------------------
        # Show thermistor depths
        # ---------------------------------------------------------
        for depth in depths:
            ax.axhline(depth,color='k',linestyle='--',linewidth=0.8,alpha=0.5)
        # ---------------------------------------------------------
        # Show periods with missing thermistor data
        # ---------------------------------------------------------
        for therm, depth in zip(therms, depths):
            missing = Df_Therm_clean[therm].isna()
            ax.scatter(dates[missing],np.full(missing.sum(), depth),color='grey',s=3,marker='s',zorder=10)
        # ---------------------------------------------------------
        # Axes
        # ---------------------------------------------------------
        ax.set_ylabel('Depth (m)', fontsize=16)
        ax.set_xlabel('Date', fontsize=16)
        ax.set_title(f'Taconnaz 2026 - Borehole {ID} - Altitude {Alt:.0f}m',fontsize=18)
        # One date tick per week
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
        plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
        # ---------------------------------------------------------
        # Colorbar
        # ---------------------------------------------------------
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('Temperature (°C)',fontsize=16)
        plt.tight_layout()

        ########################################################################################################
        ### WORK WITH MANUALLY CONSTRUCTED FILES OF STEADY TEMPERATURE FOR COMPARISON TO PREVIOUS CAMPAIGNS  ###
        ########################################################################################################
        FileName_2026_SteadyT = f'Forage_{ID.zfill(2)}.csv'
        ###Open file for past campaign that gather dat for all boreholes
        Df_2026_SteadyT = pd.read_csv(Pathroot_Steady.joinpath(FileName_2026_SteadyT), skiprows=1, header=None,names=ColNames_2026_SteadyT, sep=',', na_values=["NAN"])
        ### Force to zero slightly positive T
        Df_2026_SteadyT.loc[Df_2026_SteadyT['T'] > 0, 'T'] = 0
        ###For past campaign get the lines corresponding to Borehole of interest
        BoreholeName = f'Forage_{ID.zfill(2)}'
        Df_Past_SteadyT_ID = Df_Past_SteadyT[Df_Past_SteadyT['Borehole'].str.startswith(BoreholeName)]
        ###If requested, store basal T for each borehole
        if InterpBasalT:
            DeepestSensor2026 = Df_2026_SteadyT.loc[Df_2026_SteadyT['Sensor_Depth'].idxmin()]
            DeepestSensors2026.append(DeepestSensor2026)

        ####DO THE PLOT OF T STEADY FOR THE VARIOUS YEARS
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.set_title(f'Borehole {ID.zfill(2)} - Altitude {Alt:.0f}m', fontsize=18)

        for year, df_year in Df_Past_SteadyT_ID.groupby('Year'):
            # Sort by depth to connect points in the correct order
            df_year = df_year.sort_values('Sensor_Depth')
            ax.plot(df_year['T'],df_year['Sensor_Depth'],marker='o',markersize = 9, linewidth= 1.9, color=YearColors[year], label=str(year))
        ###Plot result for 2026
        # Sort by depth to connect points in the correct order
        Df_2026_SteadyT = Df_2026_SteadyT.sort_values('Sensor_Depth')
        ax.plot(Df_2026_SteadyT['T'], Df_2026_SteadyT['Sensor_Depth'], marker='o', markersize=9, linewidth=1.9, color=YearColors[2026],label=str(2026))


        ax.set_xlabel('Temperature (°C)', fontsize=16)
        ax.set_ylabel('Depth (m)', fontsize=16)

        ax.tick_params(labelsize=14)
        ax.grid(alpha=0.3)

        ax.legend(title='Year', fontsize=13, title_fontsize=14)

        plt.tight_layout()


    #######################################################################
    ### PLOT MAP OF INTERPOLATED BASAL TEMPERATURE FOR 2026 IF REQUIRED ###
    #######################################################################
    if InterpBasalT:
        # Build dataframe from stored basal T for 2026
        Df_DeepestSensors2026 = pd.DataFrame(DeepestSensors2026).reset_index(drop=True)
        ### Add borehole 6 which was temperated in 2014 and logically still is
        Borehole06 = {
            'Borehole': 'Forage_06',
            'Year': 2014,
            'X': 949166,
            'Y': 2104816,
            'Alt_Surf': 3649.0,
            'Borehole_Depth': 92.75,
            'Alt_Bed': 3556.25,
            'Sensor_ID': np.nan,
            'Alt_Therm': 3649.0,
            'Sensor_Depth': -92.75,
            'T': 0.0
        }
        ##Add borehole 6 in both 2026 and 2021 dataframe
        for df in [Df_DeepestSensors2026, Df_DeepestSensors2021]:
            df.loc[len(df)] = Borehole06

        ### For results of 2026, add borehole 2025 Haut which has beeen measured the summer before
        Df_DeepestSensors2026.loc[len(Df_DeepestSensors2026)] = {
            'Borehole': '25-H',
            'Year': 2025,
            'X': 949490,
            'Y': 2104877,
            'Alt_Surf': 3636.0,
            'Borehole_Depth': 85,
            'Alt_Bed': 3551,
            'Sensor_ID': np.nan,
            'Alt_Therm': 3649.0,
            'Sensor_Depth': -85,
            'T': -5.5
        }


        ###Now do same figure for 2021 and 2026 either on 2 separated figures or on a subplot
        if PlotAsSubplots:
            fig = plt.figure(figsize=(16, 8))
            gs = fig.add_gridspec(1, 3,width_ratios=[1, 1, 0.05],wspace=0.08)
            axs = [fig.add_subplot(gs[0, 0]),fig.add_subplot(gs[0, 1])]
            cax = fig.add_subplot(gs[0, 2])
        for i, (Year,Df_DeepestSensors) in enumerate([(2026, Df_DeepestSensors2026),(2021, Df_DeepestSensors2021)]):
            # Coordinates and basal temperatures
            x = Df_DeepestSensors['X'].to_numpy()
            y = Df_DeepestSensors['Y'].to_numpy()
            T = Df_DeepestSensors['T'].to_numpy()
            # Regular grid
            nx, ny = 200, 200
            x_grid = np.linspace(x.min(), x.max(), nx)
            y_grid = np.linspace(y.min(), y.max(), ny)
            X_grid, Y_grid = np.meshgrid(x_grid, y_grid)

            # Linear interpolation
            T_grid = griddata((x, y),T,(X_grid, Y_grid),method='linear')
            # ###Check the interpolation
            # print(f'\n===== {Year} =====')
            # for _, row in Df_DeepestSensors.iterrows():
            #     # Nearest grid point
            #     dist2 = (X_grid - row['X']) ** 2 + (Y_grid - row['Y']) ** 2
            #     iy, ix = np.unravel_index(np.argmin(dist2), dist2.shape)
            #
            #     print(
            #         row['Borehole'],
            #         f"Obs = {row['T']:.2f} °C",
            #         f"Grid = {T_grid[iy, ix]:.2f} °C",
            #         f"distance = {np.sqrt(dist2[iy, ix]):.1f} m"
            #     )

            ### Plot the map
            # Create/select axis
            if PlotAsSubplots:
                ax = axs[i]
            else:
                fig, ax = plt.subplots(figsize=(8, 8))
            ax.set_title(f'Obs - {Year}', fontsize=18)
            # Orthophoto
            ax.imshow(ortho_crop,extent=extent_km,origin='upper',zorder=0)
            # Temperature limits
            Tmin = -11.0
            Tmax = 0.0
            # Filled temperature contours
            levels = np.linspace(Tmin, Tmax, 12)
            cf = ax.contourf(X_grid / 1000,Y_grid / 1000,T_grid,levels=levels,cmap='RdYlBu_r',alpha=0.85,extend='min',zorder=2)
            # Temperature isolines
            contour_levels = np.arange(np.ceil(Tmin), 0.1, 1.0)
            cs = ax.contour(X_grid / 1000,Y_grid / 1000,T_grid,levels=contour_levels,colors='k',linewidths=0.8,alpha=0.8,zorder=3)
            # Mask defining the interpolated area
            interp_mask = np.isfinite(T_grid).astype(float)
            # Dashed contour around the interpolated area
            ax.contour(X_grid / 1000,Y_grid / 1000,interp_mask,levels=[0.5],colors='black',linestyles='--',linewidths=2.0,zorder=5)

            ###############################
            ### SHOW BOREHOLE LOCATIONS ###
            ###############################
            # Boreholes of currently considered year
            for _, row in Df_DeepestSensors.iterrows():
                borehole = row['Borehole']
                if borehole == '25-H':
                    ax.annotate(BoreholeLabels[borehole], xy=(row['X'] / 1000, row['Y'] / 1000), xytext=(15, 2),textcoords='offset points', fontsize=11, fontweight='bold', color='black', ha='left',va='bottom', zorder=11)
                    # dot
                    ax.scatter(row['X'] / 1000, row['Y'] / 1000, color=BoreholeColors[borehole], edgecolor='black',linewidth=1, s=240, marker='.', zorder=10)
                else:
                    # Star
                    ax.scatter(row['X'] / 1000,row['Y'] / 1000,color=BoreholeColors[borehole],edgecolor='black',linewidth=1,s=200,marker='*',zorder=10)
                    if borehole in ['Forage_03','Forage_06','Forage_10']:
                        # Borehole ID
                        ax.annotate(BoreholeLabels[borehole],xy=(row['X'] / 1000, row['Y'] / 1000),xytext=(-20, 4), textcoords='offset points',fontsize=11, fontweight='bold',color='black',ha='left',va='bottom',zorder=11)
                    else:
                                            # Borehole ID
                        ax.annotate(BoreholeLabels[borehole],xy=(row['X'] / 1000, row['Y'] / 1000),xytext=(4, 4), textcoords='offset points',fontsize=11, fontweight='bold',color='black',ha='left',va='bottom',zorder=11)
            ### Also show boreholes from the other year
            if Year == 2026:
                Df_OtherYear = Df_DeepestSensors2021
                OtherYear = 2021
            elif Year == 2021:
                Df_OtherYear = Df_DeepestSensors2026
                OtherYear = 2026

            for _, row in Df_OtherYear.iterrows():
                borehole = row['Borehole']
                # If borehole already exists in considered year, skip it,
                # EXCEPT Forage_03 because its coordinates changed (other boreholes are really close to each other between 2021 and 2026)
                if (borehole in Df_DeepestSensors['Borehole'].values and borehole != 'Forage_03'):
                    continue
                # Strong transparency for boreholes from the other year
                alpha_other = 0.42
                if borehole == '25-H':
                    ax.scatter(row['X'] / 1000,row['Y'] / 1000,color=BoreholeColors[borehole],edgecolor='black',linewidth=1,s=240,marker='.',alpha=alpha_other,zorder=10)

                else:
                    ax.scatter(row['X'] / 1000,row['Y'] / 1000,color=BoreholeColors[borehole], edgecolor='black',linewidth=1,s=200,marker='*',alpha=alpha_other,zorder=10)

                # Special annotation for the old/new position of Forage_03
                # if borehole == 'Forage_03':
                #     ax.annotate(f'03-{OtherYear}',xy=(row['X'] / 1000, row['Y'] / 1000),xytext=(-25, 4),textcoords='offset points',fontsize=10,fontweight='bold',color='black',alpha=0.5,ha='left',va='bottom',zorder=11)
            # Colorbar
            if not PlotAsSubplots:
                cbar = fig.colorbar(cf,ax=ax,pad=0.02,shrink=0.85)
                cbar.set_label(r'Basal Temperature [$^\circ$C]',fontsize=14)
                cbar.ax.tick_params(labelsize=11)
            # Map appearance
            ax.set_aspect('equal')
            ax.set_xlabel('Easting [km]', fontsize=14)
            if ((not PlotAsSubplots) or (i ==0)):
                ax.set_ylabel('Northing [km]', fontsize=14)
            else:
                ax.tick_params(axis='y', left=False, labelleft=False)

            ax.tick_params(labelsize=12)

            ax.grid(True,linewidth=0.5,alpha=0.4,zorder=1)
            if not PlotAsSubplots:
                plt.tight_layout()
            # Labels on isolines
            ax.clabel(cs,inline=True,fontsize=9,fmt='%g',manual=True)
            # ax.clabel(cs,inline=True,fontsize=9,fmt='%g')

        if PlotAsSubplots:
            cbar = fig.colorbar(cf,cax=cax,ticks=np.arange(Tmin, Tmax + 1, 1))
            cbar.set_label(r'Basal Temperature [$^\circ$C]',fontsize=14)
            cbar.ax.tick_params(labelsize=11)
            plt.tight_layout()


    plt.show()




