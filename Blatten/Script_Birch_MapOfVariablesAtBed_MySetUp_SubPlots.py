"""
@author: jbrondex

Description:
------------
This file aims at plotting figures corresponding to map of relevant variable at bed at given time
as subplots with one column per desired day

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib.ticker import FormatStrFormatter
from matplotlib.dates import date2num
import matplotlib.gridspec as gridspec
import rasterio
from rasterio.windows import from_bounds

import pandas as pd  ###To treat the csv files

from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))

###Defined Colors optimized for color-blind people:
Orange = [230 / 255, 159 / 255, 0 / 255]
SkyBlue = [86 / 255, 180 / 255, 233 / 255]
BluishGreen = [0 / 255, 158 / 255, 115 / 255]
Yellow = [240 / 255, 228 / 255, 66 / 255]
Blue = [0 / 255, 114 / 255, 178 / 255]
Vermillion = [213 / 255, 94 / 255, 0 / 255]
ReddishPurple = [204 / 255, 121 / 255, 167 / 255]
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

def Interpolate_field_slice(df, field_name, x, z): ###Returned field interpolated from mesh grid to point (x,z) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    zi = df['Z'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, zi), field, (x, z), rescale=True)
    return result

##Function below is used to sort point of contour clockwise
def sort_clockwise(list_of_xy_coords):
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x-cx, y-cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)


if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###Parameter of Simu
    C = 0.3
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 30 ##in m/d !
 #   Timestep_size = 0.01/365 ### In years !!
    Timestep_size = 0.02/365 ### In years !!
    OutputInterval = 50 ###Execution interval of the solver producing the simu output
    OutputOfPlots = [0, 1, 3, 7] #[1, 3, 7, 24] ###At which output (count) do we want to make the plot (correspond to the day of simu)
    ###Coordinate of bed point at which we want to plot the variables vs time (this is done in the other python script)
    xplot = [2630580, 2630620.9, 2630653.6]
    yplot = [1139265, 1139221.7, 1139188.3]
    colplot = ['#00FF66','#8D5FD3','#FFD42A']
    Stars = False ## Do we want to plot stars ?
    ### Size of the stars on the plot
    Marker_Size = 18 ### Size of markers to show where vertical profiles of variables vs time are plotted
    ####Plot parameters###
    ### Contours we want to plot for debris thickness
    DebrisContour_Levels = [10, 20, 30, 40, 50, 60]
    DebrisContour_LabelPositions = [(2630.73,1139.097), (2630.587,1139.089), (2630.877,1139.018), (2630.917,1138.978), (2630.736,1138.978), (2630.853,1138.930)]
    DebrisContour_Color = 'lightcoral'
    DebrisContour_Thick = 0.9
    ###Defines bounds of the plot
    xmin,xmax = 2630250, 2631250
    ymin,ymax = 1138750, 1139500
    ###Open Data corresponding to grounded mask
    # Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
    # str_Cvalue = f"{int(C * 10):02d}"
    # Filename = 'MyBirch_RS_C{}_WithD_Cont__bed.dat'.format(str_Cvalue)
    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/2_RateAndState_NoDamage_NoContact/GlacierOut/.')
    str_Cvalue = f"{int(C * 10):02d}"
    Filename = 'MyBirch_RS_C{}_WithD_NoCont__bed.dat'.format(str_Cvalue)
    Col_Names = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub',
                 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi']

    ### START OF SCRIPT###
    ###Load glacier contour
    Pathroot_Contour = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/Visu_Bed/.')
    Filename_GlacierContour = 'ContourGlacier_FromMesh3D.dat'
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
    Pathroot_DEMs= Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/1_Mesh')
    Filename_SurfGlacierDEM = 'DEMSurfPreIce.dat'
    Df_SurfGlacierDEM = pd.read_csv(Pathroot_DEMs.joinpath(Filename_SurfGlacierDEM), names=['X', 'Y', 'Z'], sep=r'\s+')
    Filename_SurfDebrisDEM = 'DEMSurfPre.dat'
    Df_SurfDebrisDEM = pd.read_csv(Pathroot_DEMs.joinpath(Filename_SurfDebrisDEM), names=['X', 'Y', 'Z'], sep=r'\s+')

    #########################################
    ####  OPEN OUTPUT OF THE SIMULATIONS:####
    #########################################
    ###load file as dataframe
    Df_wholebed_full = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Names, sep=r'\s+')
    ###Create a column with time (in years !)
    Df_wholebed_full['Time']=Timestep_size*Df_wholebed_full['Timestep']
    ###Prepare the figure : subplots with output days as column and theta, taub/N and ub as lines
    fig1, axes = plt.subplots(nrows=3, ncols=len(OutputOfPlots), figsize=(15, 12), constrained_layout=True, sharex=True, sharey=True)
    ###Prepare another figure : subplots with output days as column and taub/N and sigma_nn as lines
    fig2, axes2 = plt.subplots(nrows=2, ncols=len(OutputOfPlots), figsize=(15, 10), constrained_layout=True, sharex=True, sharey=True)
    fig2.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1, hspace=0.02, wspace=0.05)
    ### Now start loop on desired days of output
    for j,OutputOfPlot in enumerate(OutputOfPlots):
        ### Keep only lines for timestep at which we want to make the plot
        Df_wholebed=Df_wholebed_full[Df_wholebed_full['Count']==int((OutputOfPlot/(Timestep_size*OutputInterval*365))+1)].copy()
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed['un'] = Df_wholebed['u']*Df_wholebed['nx'] + Df_wholebed['v']*Df_wholebed['ny'] + Df_wholebed['w']*Df_wholebed['nz']
        ###Deduce ub from u and un
        Df_wholebed['ub'] = np.sqrt((Df_wholebed['u']-Df_wholebed['un']*Df_wholebed['nx'])**2+(Df_wholebed['v']-Df_wholebed['un']*Df_wholebed['ny'])**2 +(Df_wholebed['w']-Df_wholebed['un']*Df_wholebed['nz'])**2 )
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed['Sigma_nn'] = (Df_wholebed['fx']*Df_wholebed['nx'] + Df_wholebed['fy']*Df_wholebed['ny'] + Df_wholebed['fz']*Df_wholebed['nz'])/Df_wholebed['nodearea']
        ###Create refined regular grid and interpolate field over grid (only once)
        if j ==0:
            xbed = np.arange(np.floor(np.min(Df_wholebed['X']))-10,np.floor(np.max(Df_wholebed['X']))+10,0.5)
            ybed = np.arange(np.floor(np.min(Df_wholebed['Y']))-10,np.floor(np.max(Df_wholebed['Y']))+10,0.5)
            Xbed, Ybed = np.meshgrid(xbed, ybed)
            ### Interpolate DEMs over same grid (just needed once)
            ZSurfGlacier = Interpolate_field(Df_SurfGlacierDEM, 'Z', Xbed, Ybed)
            ZSurfDebris =  Interpolate_field(Df_SurfDebrisDEM, 'Z', Xbed, Ybed)
            Hdebris_notclipped = ZSurfDebris - ZSurfGlacier
            ### If H debris is less than 1.5m -> Zero
            Hdebris = np.where(Hdebris_notclipped < 1.5, 0, Hdebris_notclipped)
            Hdebris_smoothed = gaussian_filter(Hdebris, sigma=6)
        taub = Interpolate_field(Df_wholebed,'taub',Xbed,Ybed)
        N = Interpolate_field(Df_wholebed,'N',Xbed,Ybed)
        ub = Interpolate_field(Df_wholebed,'ub',Xbed,Ybed)
        Sigma_nn = Interpolate_field(Df_wholebed,'Sigma_nn',Xbed,Ybed)
        theta = Interpolate_field(Df_wholebed,'theta',Xbed,Ybed)
        ##############################################################
        ###     PLOT MAPS OF VARIABLE OF INTEREST FOR FIGURE 1     ###
        ##############################################################
        ###~~~~~~~Row 1 : Theta ~~~~~~~~~~~~~~~
        ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
        ###Get the proper subplot
        ax = axes[0,j]
        if j == 0:
            ax.set_ylabel(r'Y [km]', fontsize=18)
            # ax.text(-0.15, 0.5, r'Cavitation', va='center', ha='center', rotation='vertical',fontsize=18, transform=ax.transAxes)
        ax.tick_params(labelsize=14)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('Day {}'.format(OutputOfPlot), fontsize=21, weight='bold')
        #shading
        clevs=np.arange(0.0, 1.01, 0.1)   ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'cividis'
        cmap = cmx.get_cmap(Colormap)
        #colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed/1000, Ybed/1000, 1-theta, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
        #CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, theta, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax, alpha=0.6, zorder=1)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed/1000, Ybed/1000, Hdebris_smoothed , levels= DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-')
        ###Put labels with position determined manually
        ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d", manual=DebrisContour_LabelPositions, fontsize=7.5)
        ax.plot(xc/1000, yc/1000 , color='k', linewidth=2)
        # ###Plot stars where we plot variable versus time
        if Stars:
            for xstar,ystar,colstar in zip(xplot,yplot,colplot):
                ax.plot(xstar/1000, ystar/1000, color=colstar, marker='*', markersize=Marker_Size, markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1, CS1_debriscontour]:
            c.set_clip_path(patch)
        if j == len(OutputOfPlots)-1:
            cbar = fig1.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
            cbar.set_ticks(np.arange(0, 1.01, 0.2))
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(r'$\theta$', fontsize=28)

        ###~~~~~~~Row 2 : taub/N ~~~~~~~~~~~~~~~
        ### Make the division by being careful of places where N might be zero or NaN
        taub_over_N_raw = np.divide(taub, N, out=np.full_like(taub, np.nan), where=(N != 0))
        # Mask NaNs and infs for safe plotting
        taub_over_N = np.ma.masked_invalid(taub_over_N_raw)
        ###Get the proper subplot
        ax = axes[1,j]
        if j == 0:
            ax.set_ylabel(r'Y [km]', fontsize=18)
            # ax.text(-0.15, 0.5, r'$\tau_\mathrm{b}/N$ ratio', va='center', ha='center',rotation='vertical', fontsize=18, transform=ax.transAxes)
        ax.tick_params(labelsize=14)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        #shading
        clevs=np.arange(0.0, 0.501, 0.05)   ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'YlGnBu_r'
        cmap = cmx.get_cmap(Colormap)
        #colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed/1000, Ybed/1000, taub_over_N, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
        #CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, taub_over_N,  cmap=cmap, shading='auto', vmin=vmin, vmax =vmax ,alpha=0.6,zorder=1)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed/1000, Ybed/1000, Hdebris_smoothed , levels= DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-')
        ax.plot(xc/1000, yc/1000 , color='k', linewidth=2)
        #### For the steady state of cavitation (so only if day of output is greater than day at which steady state is reached), hatched places
        #### where taub/N has reached the threshold C (with a margin of 5%)
        if j == len(OutputOfPlots)-1:
            CS1_contour = ax.contour(Xbed / 1000, Ybed / 1000, taub_over_N, levels=[C-0.05*C], colors='black', linewidths=1.2, linestyles='-')
            # # Hatch the area where taub/N => C
            CS1_hatch = ax.contourf(Xbed / 1000, Ybed / 1000, taub_over_N,levels=[C-0.05*C, np.max(taub_over_N)],colors='none', hatches=['////'], alpha=0)
        # ###Plot stars where we plot variable versus time
        if Stars:
            for xstar,ystar,colstar in zip(xplot,yplot,colplot):
                ax.plot(xstar/1000, ystar/1000, color=colstar, marker='*', markersize=Marker_Size, markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1,CS1_debriscontour]:
            c.set_clip_path(patch)
        if j == len(OutputOfPlots)-1:
            for c in [CS1_contour]:
                c.set_clip_path(patch)
            for c in [CS1_hatch]:
                c.set_clip_path(patch)
        if j == len(OutputOfPlots) - 1:
            cbar = fig1.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
            cbar.set_ticks(np.arange(0, 0.501, 0.1))
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(r'$\tau_\mathrm{b}/N$', fontsize=28)

        ###~~~~~~~Row 3 : ub ~~~~~~~~~~~~~~~
        ###Get the proper subplot
        ax = axes[2,j]
        ax.set_xlabel(r'X [km]', fontsize=18)
        if j == 0:
            ax.set_ylabel(r'Y [km]', fontsize=18)
            # ax.text(-0.15, 0.5, r'Basal velocity', va='center', ha='center', rotation='vertical', fontsize=18, transform=ax.transAxes)
        ax.tick_params(labelsize=14)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        #shading
        clevs=np.arange(0.0, velocitymax+0.1, 2.5)   ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'viridis'
        cmap = cmx.get_cmap(Colormap)
        #colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed/1000, Ybed/1000, ub/365, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
        #CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, ub/365,  cmap=cmap, shading='auto', vmin=vmin, vmax =vmax ,alpha=0.6,zorder=1)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed/1000, Ybed/1000, Hdebris_smoothed , levels= DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-')
        ax.plot(xc/1000, yc/1000 , color='k', linewidth=2)
        ###Plot stars where we plot variable versus time
        if Stars:
            for xstar,ystar,colstar in zip(xplot,yplot,colplot):
                ax.plot(xstar/1000, ystar/1000, color=colstar, marker='*', markersize=Marker_Size, markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1, CS1_debriscontour]:
            c.set_clip_path(patch)
        if j == len(OutputOfPlots) - 1:
            cbar = fig1.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
            cbar.set_ticks(np.arange(0, velocitymax+0.1, 5))
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(r'$\mathrm{u_b}$ [m/d]', fontsize=28)

    ##############################################################
    ###     PLOT MAPS OF VARIABLE OF INTEREST FOR FIGURE 2     ###
    ##############################################################
        ###~~~~~~~Row 1 : taub/N ~~~~~~~~~~~~~~~
        ###Get the proper subplot
        ax = axes2[0, j]
        if j == 0:
            ax.set_ylabel(r'Y [km]', fontsize=18)
            # ax.text(-0.15, 0.5, r'$\tau_\mathrm{b}/N$ ratio', va='center', ha='center',rotation='vertical', fontsize=18, transform=ax.transAxes)
        ax.tick_params(labelsize=14)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('Day {}'.format(OutputOfPlot), fontsize=21, weight='bold')
        # shading
        clevs = np.arange(0.0, 0.501, 0.05)  ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'YlGnBu_r'
        cmap = cmx.get_cmap(Colormap)
        # colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, taub_over_N, levels=clevs, cmap=cmap, extend='both', alpha=0.7,
                          antialiased=True)
        # CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, taub_over_N,  cmap=cmap, shading='auto', vmin=vmin, vmax =vmax ,alpha=0.6,zorder=1)
        # ### Plot contours of debris thickness
        # CS1_debriscontour = ax.contour(Xbed/1000, Ybed/1000, Hdebris_smoothed , levels= DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-')
        # ###Put labels with position determined manually
        # ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d", manual=DebrisContour_LabelPositions, fontsize=7.5)
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        #### For the steady state of cavitation (so only if day of output is greater than day at which steady state is reached), hatched places
        #### where taub/N has reached the threshold C (with a margin of 5%)
        if j == len(OutputOfPlots) - 1:
            CS1_contour = ax.contour(Xbed / 1000, Ybed / 1000, taub_over_N, levels=[C - 0.05 * C], colors='black',
                                     linewidths=1.2, linestyles='-')
            # # Hatch the area where taub/N => C
            CS1_hatch = ax.contourf(Xbed / 1000, Ybed / 1000, taub_over_N,
                                    levels=[C - 0.05 * C, np.max(taub_over_N)], colors='none', hatches=['////'],
                                    alpha=0)
        # ###Plot stars where we plot variable versus time
        if Stars:
            for xstar, ystar, colstar in zip(xplot, yplot, colplot):
                ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Marker_Size,
                        markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1]:
            c.set_clip_path(patch)
        if j == len(OutputOfPlots) - 1:
            for c in [CS1_contour]:
                c.set_clip_path(patch)
            for c in [CS1_hatch]:
                c.set_clip_path(patch)
        if j == len(OutputOfPlots) - 1:
            cbar = fig1.colorbar(CS1, ax=ax, shrink=0.6, pad=0.05)
            cbar.set_ticks(np.arange(0, 0.501, 0.1))
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(r'$\tau_\mathrm{b}/N$', fontsize=28)

        ###~~~~~~~Row 2 : Sigma_nn ~~~~~~~~~~~~~~~
        ###Get the proper subplot
        ax = axes2[1, j]
        ax.set_xlabel(r'X [km]', fontsize=18)
        if j == 0:
            ax.set_ylabel(r'Y [km]', fontsize=18)
        ax.tick_params(labelsize=14)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        # shading
        clevs = np.arange(-2.0, 2.01, 0.05)  ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'RdBu'
        cmap = cmx.get_cmap(Colormap)
        # colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        # CS1 = ax.contourf(Xbed/1000, Ybed/1000, Sigma_nn, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
        CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, Sigma_nn, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax, alpha=0.6,
                            zorder=1)
        # ### Plot contours of debris thickness
        # CS1_debriscontour = ax.contour(Xbed/1000, Ybed/1000, Hdebris_smoothed , levels= DebrisContour_Levels, colors=DebrisContour_Color, linewidths=DebrisContour_Thick, linestyles='-')
        # ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)

        #### For the steady state of cavitation (so only if day of output is greater than day at which steady state is reached), hatched places
        #### where taub/N has reached the threshold C (with a margin of 5%)
        if j == len(OutputOfPlots) - 1:
            CS1_contour = ax.contour(Xbed / 1000, Ybed / 1000, taub_over_N, levels=[C - 0.05 * C], colors='black', linewidths=1.2, linestyles='-')
        # ###Plot stars where we plot variable versus time
        if Stars:
            for xstar, ystar, colstar in zip(xplot, yplot, colplot):
                ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Marker_Size, markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1]:
            c.set_clip_path(patch)
        if j == len(OutputOfPlots) - 1:
            for c in [CS1_contour]:
                c.set_clip_path(patch)
        if j == len(OutputOfPlots) - 1:
            cbar = fig1.colorbar(CS1, ax=ax, shrink=0.6, pad=0.05)
            cbar.set_ticks(np.arange(-2, 2.01, 0.5))
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=28)


    plt.show()




