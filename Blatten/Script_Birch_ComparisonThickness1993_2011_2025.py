"""
@author: jbrondex

Description:
------------
This file aims at plotting figures comparing the thickness and contour of Birch over several years

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.lines import lineStyles
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch
from shapely.geometry import Polygon
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


# Function to add label along a contour
def label_along_line(ax, x, y, text, idx=None, offset=0.02, **kwargs):
    """
    Places text along a line at index `idx` (default: middle),
    rotated to follow the local tangent direction.
    """
    if idx is None:
        idx = len(x) // 2

    # Compute tangent angle
    dx = x[idx+1] - x[idx-1]
    dy = y[idx+1] - y[idx-1]
    angle = np.degrees(np.arctan2(dy, dx))
    # Normalize tangent and compute normal
    tangent = np.array([dx, dy])
    tangent /= np.linalg.norm(tangent)
    normal = np.array([-tangent[1], tangent[0]])  # Rotate tangent by +90°
    # Offset point along normal
    x_text = x[idx] + offset * normal[0]
    y_text = y[idx] + offset * normal[1]

    ax.text(x_text, y_text, text,
            rotation=angle,
            rotation_mode='anchor',
            ha='center', va='center',
            **kwargs)
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
    ###Parameter of Script
    Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/1_Mesh/.')
    Years = [1993,2011,2025]
    Col_Name = ['Count', 'BC', 'X', 'Y', 'Z', 'H']
    ####Plot parameters###
    ###Defines bounds of the plot
    xmin,xmax = 2630250, 2631250
    ymin,ymax = 1138750, 1139500
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

    ### START OF SCRIPT###
    xbed = np.arange(xmin + 30, xmax - 30, 0.5)
    ybed = np.arange(ymin + 30, ymax - 30, 0.5)
    Xbed, Ybed = np.meshgrid(xbed, ybed)
    ###Prepare the figure : subplots with absolute thickness for each year
    fig1, axes1 = plt.subplots(nrows=1, ncols=len(Years), figsize=(12, 8), constrained_layout=True, sharex=True, sharey=True)
    ###Start loop on years
    Allcontours=pd.DataFrame()
    Allthickness=pd.DataFrame()
    for i,Year in enumerate(Years):
        ###Load glacier contour
        if Year == 2025:
            Filename_GlacierContour = 'ContourGlacier_FromMesh3D.dat'
            Pathroot_Contour = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/Visu_Bed/.')
        else:
            Filename_GlacierContour = 'Contour_Birch_{}.dat'.format(Year)
            Pathroot_Contour = Pathroot.joinpath('MESH{}'.format(Year))
        ###Print some info
        print('Opening file:', Pathroot_Contour.joinpath(Filename_GlacierContour))
        ###Load glacier contour
        Df_GlacierContour = pd.read_csv(Pathroot_Contour.joinpath(Filename_GlacierContour), names=['X', 'Y'], sep=r'\s+')
        # Close the contour: Copy the first row and append it to the end
        Df_GlacierContour = pd.concat([Df_GlacierContour, Df_GlacierContour.iloc[[0]]], ignore_index=True)
        ##Add a column corresponding to year
        Df_GlacierContour['Year'] = Year
        ##Store x and y of contour of current year
        xc = Df_GlacierContour['X'].values
        yc = Df_GlacierContour['Y'].values
        ###Fill the Allcontours dataframe
        Allcontours = pd.concat([Allcontours, Df_GlacierContour], ignore_index=True)
        ###Load Ice thickness
        if Year == 2025:
            Filename_Thick = './MESH{}/GlacierOut/MyBirch_Init_DEMCorrected_GL_NoFS__bed.dat'.format(Year)
        else:
            Filename_Thick = './MESH{}/GlacierOut/MyBirch{}_Init_DEMCorrected_GL_NoFS__bed.dat'.format(Year,Year)
        Df_Thick = pd.read_csv(Pathroot.joinpath(Filename_Thick), names=Col_Name, sep=r'\s+')
        ##Add a column corresponding to year
        Df_Thick['Year'] = Year
        ###Fill the Allthickness dataframe
        Allthickness = pd.concat([Allthickness, Df_Thick], ignore_index=True)
        ###Interpolate thickness over Xbed Ybed:
        ThickGlacier = Interpolate_field(Df_Thick, 'H', Xbed, Ybed)
        ####Replace NaN by zeros : No Ice
        ThickGlacier[np.isnan(ThickGlacier)] = 0.0
        ##############################################################
        ###     PLOT MAPS OF VARIABLE OF INTEREST FOR FIGURE 1     ###
        ##############################################################
        ###~~~~~~~Row 1 : Theta ~~~~~~~~~~~~~~~
        ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
        ###Get the proper subplot
        ax = axes1[i]
        ax.set_xlabel(r'X [km]', fontsize=18)
        if i== 0:
            ax.set_ylabel(r'Y [km]', fontsize=18)
            # ax.text(-0.15, 0.5, r'Cavitation', va='center', ha='center', rotation='vertical',fontsize=18, transform=ax.transAxes)
        ax.tick_params(labelsize=14)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('{}'.format(Year), fontsize=21, weight='bold')
        #shading
        clevs=np.arange(0.0, 60, 2)   ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'inferno'
        cmap = cmx.get_cmap(Colormap)
        #colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed/1000, Ybed/1000, ThickGlacier, levels=clevs, cmap=cmap, extend='both', alpha=0.9, antialiased=True)
        ax.plot(xc/1000, yc/1000 , color='k', linewidth=2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1]:
            c.set_clip_path(patch)
        if i == len(Years)-1:
            cbar = fig1.colorbar(CS1, ax=ax, shrink=0.45, pad=0.05)
            cbar.set_ticks(np.arange(0, 60.01, 10))
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(r'H[m]', fontsize=28)

    ######################################
    ###     PLOT DH/DT IN FIGURE 2     ###
    ######################################
    ###Prepare the figure : subplots with absolute thickness for each year
    fig2, axes2 = plt.subplots(nrows=1, ncols=len(Years), figsize=(12, 8), constrained_layout=True, sharex=True, sharey=True)
    ###Get thickness on regular grid year per year
    ThickGlacier2025 = Interpolate_field(Allthickness[Allthickness['Year']==2025], 'H', Xbed, Ybed)
    ThickGlacier2011 = Interpolate_field(Allthickness[Allthickness['Year']==2011], 'H', Xbed, Ybed)
    ThickGlacier1993 = Interpolate_field(Allthickness[Allthickness['Year']==1993], 'H', Xbed, Ybed)
    ####Replace NaN by zeros : No Ice
    ThickGlacier2025[np.isnan(ThickGlacier2025)] = 0.0
    ThickGlacier2011[np.isnan(ThickGlacier2011)] = 0.0
    ThickGlacier1993[np.isnan(ThickGlacier1993)] = 0.0
    ###Get contour year per year
    Contour2025 = Allcontours[Allcontours['Year']==2025]
    Contour2011 = Allcontours[Allcontours['Year']==2011]
    Contour1993 = Allcontours[Allcontours['Year']==1993]

    ###Subplot 1 is 2011-1993
    ax=axes2[0]
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Y [km]', fontsize=18)
    ax.tick_params(labelsize=14)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('dh 2011-1993', fontsize=21, weight='bold')
    # shading
    clevs = np.arange(-30, 30, 1)  ## cbar for shading
    vmin, vmax = clevs[0], clevs[-1]
    Colormap = 'coolwarm_r'
    cmap = cmx.get_cmap(Colormap)
    # colorbar
    # levs_ticks = [1, 10, 100, 1000, 10000]
    # --- Plot orthophoto as background
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
    ###Fills up the map with colors for SigmaEq
    CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, ThickGlacier2011-ThickGlacier1993, levels=clevs, cmap=cmap, extend='both', alpha=0.9,antialiased=True)
    ax.plot(Contour1993['X'].values / 1000, Contour1993['Y'].values / 1000, color='k', linewidth=2, linestyle = '--')
    ax.plot(Contour2011['X'].values / 1000, Contour2011['Y'].values / 1000, color='k', linewidth=2, linestyle = ':')
    # Add labels along contours
    label_along_line(ax, Contour2011['X'].values / 1000, Contour2011['Y'].values / 1000, '2011', idx=len(Contour2011['X'])-5, fontsize=12, color='k')
    label_along_line(ax, Contour1993['X'].values / 1000, Contour1993['Y'].values / 1000, '1993', idx=len(Contour1993['X'])-12, fontsize=12, color='k')

    ###Below we remove colors that are outside of glacier contour
    clippath = mpltPath(np.c_[ Contour1993['X'].values/ 1000, Contour1993['Y'].values / 1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS1]:
        c.set_clip_path(patch)

    ###Subplot 2 is 2025-2011
    ax=axes2[1]
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.tick_params(labelsize=14)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('dh 2025-2011', fontsize=21, weight='bold')
    # shading
    clevs = np.arange(-30, 30, 1)  ## cbar for shading
    vmin, vmax = clevs[0], clevs[-1]
    Colormap = 'coolwarm_r'
    cmap = cmx.get_cmap(Colormap)
    # colorbar
    # levs_ticks = [1, 10, 100, 1000, 10000]
    # --- Plot orthophoto as background
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
    ###Fills up the map with colors for SigmaEq
    CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, ThickGlacier2025-ThickGlacier2011, levels=clevs, cmap=cmap, extend='both', alpha=0.9,
                      antialiased=True)
    ax.plot(Contour2025['X'].values / 1000, Contour2025['Y'].values / 1000, color='k', linewidth=2)
    ax.plot(Contour2011['X'].values / 1000, Contour2011['Y'].values / 1000, color='k', linewidth=2, linestyle = '--')
    # Add labels along contours
    label_along_line(ax, Contour2011['X'].values / 1000, Contour2011['Y'].values / 1000, '2011', idx=len(Contour2011['X'])-5, fontsize=12, color='k')
    label_along_line(ax, Contour2025['X'].values / 1000, Contour2025['Y'].values / 1000, '2025', idx=len(Contour2025['X'])-140, fontsize=12, color='k')
    ###Below we remove colors that are outside of glacier contour
    clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS1]:
        c.set_clip_path(patch)

    ###Subplot 3 is 2025-1993
    ax = axes2[2]
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Y [km]', fontsize=18)
    ax.tick_params(labelsize=14)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('dh 2025-1993', fontsize=21, weight='bold')
    # shading
    clevs = np.arange(-30, 30, 1)  ## cbar for shading
    vmin, vmax = clevs[0], clevs[-1]
    Colormap = 'coolwarm_r'
    cmap = cmx.get_cmap(Colormap)
    # colorbar
    # levs_ticks = [1, 10, 100, 1000, 10000]
    # --- Plot orthophoto as background
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
    ###Fills up the map with colors for SigmaEq
    CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, ThickGlacier2025 - ThickGlacier1993, levels=clevs, cmap=cmap,
                      extend='both', alpha=0.9,
                      antialiased=True)
    ax.plot(Contour2025['X'].values / 1000, Contour2025['Y'].values / 1000, color='k', linewidth=2)
    ax.plot(Contour1993['X'].values / 1000, Contour1993['Y'].values / 1000, color='k', linewidth=2, linestyle='--')
    # Add labels along contours
    label_along_line(ax, Contour2025['X'].values / 1000, Contour2025['Y'].values / 1000, '2025', idx=len(Contour2025['X'])-140, fontsize=12, color='k')
    label_along_line(ax, Contour1993['X'].values / 1000, Contour1993['Y'].values / 1000, '1993', idx=len(Contour1993['X'])-10, offset=-0.03, fontsize=12, color='k')
    ###Below we remove colors that are outside of the wider of the glacier contour
    clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS1]:
        c.set_clip_path(patch)
    cbar = fig2.colorbar(CS1, ax=ax, shrink=0.45, pad=0.05)
    cbar.set_ticks(np.arange(-30, 30.01, 10))
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label(r'dH[m]', fontsize=28)


    ###Make a standalone fig for dh 2025-1993 with adapted colorbar scale
    # === Create a standalone figure and axis ===
    fig, ax = plt.subplots(figsize=(8, 7))  # choose suitable size

    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Y [km]', fontsize=18)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.5)

    # Force same aspect ratio
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('dh 2025-1993', fontsize=21, weight='bold')

    # --- Color levels and colormap ---
    clevs = np.arange(-45, 45, 1)
    vmin, vmax = clevs[0], clevs[-1]
    Colormap = 'coolwarm_r'
    cmap = cmx.get_cmap(Colormap)

    # --- Background orthophoto ---
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)

    # --- Filled contour (dh) ---
    CS1 = ax.contourf(Xbed / 1000, Ybed / 1000,
                      ThickGlacier2025 - ThickGlacier1993,
                      levels=clevs, cmap=cmap,
                      extend='both', alpha=0.9, antialiased=True)

    # --- Glacier outlines ---
    ax.plot(Contour2025['X'].values / 1000, Contour2025['Y'].values / 1000,
            color='k', linewidth=2)
    ax.plot(Contour1993['X'].values / 1000, Contour1993['Y'].values / 1000,
            color='k', linewidth=2, linestyle='--')

    # --- Labels along contours ---
    label_along_line(ax, Contour2025['X'].values / 1000, Contour2025['Y'].values / 1000,
                     '2025', idx=len(Contour2025['X']) - 140,
                     fontsize=12, color='k')
    label_along_line(ax, Contour1993['X'].values / 1000, Contour1993['Y'].values / 1000,
                     '1993', idx=len(Contour1993['X']) - 10,
                     offset=-0.03, fontsize=12, color='k')

    # --- Clip to glacier boundary ---
    clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    CS1.set_clip_path(patch)

    # --- Colorbar ---
    cbar = fig.colorbar(CS1, ax=ax, shrink=0.8, pad=0.05)
    cbar.set_ticks(np.arange(-45, 45.01, 15))
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label(r'dH [m]', fontsize=28)

    plt.tight_layout()

    plt.show()




