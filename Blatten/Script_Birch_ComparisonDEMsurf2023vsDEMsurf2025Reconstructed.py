"""
@author: jbrondex

Description:
------------
This file aims at plotting figures comparing the surface DEM reconstructed by Adrien combining the 23rd of May 25 and
the 2023 surface DEM and the pre-event surface DEM reconstructed by Mylène Jacquemart assuming mass conservation from
2023 and dynamic transfer from back to front. I also show the debris thickness produced in each cases.

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
from scipy.spatial.distance import cdist
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
    ####Plot parameters###
    ###Defines bounds of the plot
    xmin,xmax = 2630250, 2631250
    ymin,ymax = 1138750, 1139500

    ### Contours we want to plot for debris thickness
    DebrisContour_Levels = [10, 20, 30, 40, 50, 60]
    DebrisContour_LabelPositions_Adrien = [(2630.73,1139.097), (2630.587,1139.089), (2630.877,1139.018), (2630.917,1138.978), (2630.736,1138.978), (2630.853,1138.930)]
    DebrisContour_LabelPositions_Mylene = [(2630.801,1139.1545), (2630.78,1139.079), (2630.899,1139.043), (2630.8928,1139.006), (2630.899,1138.95), (2630.801,1138.97)]
    DebrisContour_Thick = 0.9
    ### Not Same colour for debris contour in olver vs new approach (i.e. using Adrien's DEM versus Mylene's DEM for glacier surface)
    DebrisContour_Color_Adrien = 'r' #'grey'
    DebrisContour_Color_Mylene = 'k'

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

    ###Load glacier contour
    Pathroot_Contour = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Avr2026/1_Mesh/MESH2025/.')
    Filename_GlacierContour = 'Contour_Birch_2025_Julien.dat'
    ###Load glacier contour
    Df_GlacierContour = pd.read_csv(Pathroot_Contour.joinpath(Filename_GlacierContour), names=['X', 'Y'], sep=r'\s+')
    # Close the contour: Copy the first row and append it to the end
    Df_GlacierContour = pd.concat([Df_GlacierContour, Df_GlacierContour.iloc[[0]]], ignore_index=True)
    ##Store x and y of contour of current year
    xc = Df_GlacierContour['X'].values
    yc = Df_GlacierContour['Y'].values

    ### Load DEM of surface of glacier and rock deposit for Adrien and Mylene reconstructed DEMs
    # Pathroot_Adrien = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_Jan2026/1_Mesh/MESH2025')
    # Filename_SurfGlacierDEM_Adrien = 'DEM_surf_Corrected.dat'
    Pathroot_Adrien = Path('/home/brondexj/Documents/Expertises/Blatten/DEMs/.')
    Filename_SurfGlacierDEM_Adrien = 'DEM_surf_2024.dat'
    Df_SurfGlacierDEM_Adrien = pd.read_csv(Pathroot_Adrien.joinpath(Filename_SurfGlacierDEM_Adrien), names=['X', 'Y', 'Z'], sep=r'\s+')

    Pathroot_Mylene = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Avr2026/1_Mesh/MESH2025')
    Filename_SurfGlacierDEM_Mylene = 'DEM_surf_2025_Reconstructed_FromAdrien_MJacquemart_Smoothed.dat'
    Df_SurfGlacierDEM_Mylene = pd.read_csv(Pathroot_Mylene.joinpath(Filename_SurfGlacierDEM_Mylene), names=['X', 'Y', 'Z'], sep=r'\s+')

    Filename_SurfDebrisDEM = 'DEMSurfPre.dat'
    Df_SurfDebrisDEM = pd.read_csv(Pathroot_Mylene.joinpath(Filename_SurfDebrisDEM), names=['X', 'Y', 'Z'], sep=r'\s+')

    ### START OF SCRIPT###
    xbed = np.arange(xmin + 30, xmax - 30, 0.5)
    ybed = np.arange(ymin + 30, ymax - 30, 0.5)
    Xbed, Ybed = np.meshgrid(xbed, ybed)
    ### Interpolate DEMs over same grid for debris thickness distribution (just needed once)
    ZSurfGlacier_Adrien = Interpolate_field(Df_SurfGlacierDEM_Adrien, 'Z', Xbed, Ybed)
    ZSurfGlacier_Mylene = Interpolate_field(Df_SurfGlacierDEM_Mylene, 'Z', Xbed, Ybed)
    ZSurfDebris = Interpolate_field(Df_SurfDebrisDEM, 'Z', Xbed, Ybed)

    Hdebris_notclipped_Adrien = ZSurfDebris - ZSurfGlacier_Adrien
    ### For contour lines : If H debris is less than 1.0m -> Zero
    Hdebris_Adrien = np.where(Hdebris_notclipped_Adrien < 1.0, 0, Hdebris_notclipped_Adrien)
    Hdebris_smoothed_Adrien  = gaussian_filter(Hdebris_Adrien , sigma=6)

    Hdebris_notclipped_Mylene = ZSurfDebris - ZSurfGlacier_Mylene
    ### For contour lines : If H debris is less than 1.0m -> Zero
    Hdebris_Mylene = np.where(Hdebris_notclipped_Mylene < 1.0, 0, Hdebris_notclipped_Mylene)
    Hdebris_smoothed_Mylene  = gaussian_filter(Hdebris_Mylene , sigma=6)

    #########FIG 1 : DEBRIS THICKNESS IN EACH APPROACH ############
    ###Prepare the figure : subplots with absolute debris thickness for the 2 approaches
    fig1, axes1 = plt.subplots(nrows=1, ncols=2, figsize=(12, 8), constrained_layout=True, sharex=True, sharey=True)
    ###Get the proper subplot
    ax = axes1[0]
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Y [km]', fontsize=18)
    ax.tick_params(labelsize=14)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect('equal', adjustable='box')
    # ax.set_title('Surface DEM of Adrien', fontsize=21, weight='bold')
    ax.set_title('Surface 2024', fontsize=21, weight='bold')
    #shading
    clevs=np.arange(0.0, 70.01, 5)   ## cbar for shading
    vmin, vmax = clevs[0], clevs[-1]
    Colormap = 'inferno' #'cividis'
    cmap = cmx.get_cmap(Colormap)
    #colorbar
    # levs_ticks = [1, 10, 100, 1000, 10000]
    # --- Plot orthophoto as background
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
    ###Fills up the map with colors for debris thickness
    CS1 = ax.contourf(Xbed/1000, Ybed/1000, Hdebris_notclipped_Adrien, levels=clevs, cmap=cmap, extend='both', alpha=0.9, antialiased=True)
    ax.plot(xc/1000, yc/1000 , color='k', linewidth=2)
    ###Below we remove colors that are outside of glacier contour
    clippath = mpltPath(np.c_[xc/1000, yc/1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS1]:
        c.set_clip_path(patch)

    ax = axes1[1]
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.tick_params(labelsize=14)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect('equal', adjustable='box')
    # ax.set_title('Surface DEM of Mylene', fontsize=21, weight='bold')
    ax.set_title('Surface Currently Used', fontsize=21, weight='bold')
    # shading
    clevs = np.arange(0.0, 70.01, 5)  ## cbar for shading
    vmin, vmax = clevs[0], clevs[-1]
    # colorbar
    # levs_ticks = [1, 10, 100, 1000, 10000]
    # --- Plot orthophoto as background
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
    ###Fills up the map with colors for debris thickness
    CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, Hdebris_notclipped_Mylene, levels=clevs, cmap=cmap, extend='both',alpha=0.9, antialiased=True)
    ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    ###Below we remove colors that are outside of glacier contour
    clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS1]:
        c.set_clip_path(patch)
    cbar = fig1.colorbar(CS1, ax=ax, shrink=0.45, pad=0.05)
    cbar.set_ticks(np.arange(0, 70.01, 10))
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label(r'Hdebris[m]', fontsize=28)



    ###Make a standalone fig for dh between two approaches with contour lines of debris
    # === Create a standalone figure and axis ===
    fig, ax = plt.subplots(figsize=(8, 7))  # choose suitable size

    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Y [km]', fontsize=18)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.5)

    # Force same aspect ratio
    ax.set_aspect('equal', adjustable='box')
    # ax.set_title('Surface Mylene - Surface Adrien', fontsize=21, weight='bold')
    ax.set_title('Surface Currently Used - Surface 2024', fontsize=21, weight='bold')
    # --- Color levels and colormap ---
    clevs = np.arange(-25, 25, 1)
    vmin, vmax = clevs[0], clevs[-1]
    Colormap = 'coolwarm_r'
    cmap = cmx.get_cmap(Colormap)

    # --- Background orthophoto ---
    ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)

    # --- Filled contour (dh) ---
    CS1 = ax.contourf(Xbed / 1000, Ybed / 1000,ZSurfGlacier_Mylene - ZSurfGlacier_Adrien,levels=clevs, cmap=cmap,extend='both', alpha=0.9, antialiased=True)
    # --- Glacier outlines ---
    ax.plot(xc / 1000, yc / 1000,color='k', linewidth=2)
    ### Plot contours of debris thickness
    CS1_debriscontour_Adrien = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed_Adrien, levels=DebrisContour_Levels,colors=DebrisContour_Color_Adrien, linewidths=DebrisContour_Thick, linestyles='-')
    CS1_debriscontour_Mylene = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed_Mylene,levels=DebrisContour_Levels, colors=DebrisContour_Color_Mylene, linewidths=DebrisContour_Thick, linestyles=':')
    ###Put labels with position determined manually
    ax.clabel(CS1_debriscontour_Adrien, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d",manual=DebrisContour_LabelPositions_Adrien, fontsize=7.5)
    ax.clabel(CS1_debriscontour_Mylene, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d",manual=DebrisContour_LabelPositions_Mylene, fontsize=7.5)
    # --- Clip to glacier boundary ---
    clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS1, CS1_debriscontour_Adrien, CS1_debriscontour_Mylene]:
        c.set_clip_path(patch)

    # --- Colorbar ---
    cbar = fig.colorbar(CS1, ax=ax, shrink=0.8, pad=0.05)
    cbar.set_ticks(np.arange(-25, 25.01, 5))
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label(r'd$z_\mathrm{s}$ [m]', fontsize=28)

    plt.tight_layout()

    plt.show()




