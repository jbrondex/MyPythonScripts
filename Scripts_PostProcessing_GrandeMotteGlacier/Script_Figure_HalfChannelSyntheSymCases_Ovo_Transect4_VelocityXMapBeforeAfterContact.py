"""
@author: jbrondex

Description:
------------
This script is used to produce figure showing a map of velocity X around channel (synthetic sym) only for Transect 4 and shape ovoide right before and right after
contact
"""


################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.patches as Patch

from pathlib import Path

import numpy as np
from matplotlib.ticker import FormatStrFormatter

import pandas as pd  ###To treat the csv files
from scipy import interpolate
from scipy.interpolate import griddata

from matplotlib import ticker
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch

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
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result
################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=24)

###Fig 1 is subplots for right before and right after contact
fig, axes = plt.subplots(1, 2, figsize=(25, 50), sharey=True)
fig.tight_layout(pad=28.0, w_pad=-0.1)
# Adjust layout to make room for the legend
for j in range(2):
    ax = axes[j]
    ax.set_xlim([-1.0, 26])
    ax.set_ylim([2774, 2829.5])##en altitude
    # ax.locator_params(axis = 'y', nbins = 5)
    ####Plot vertical lines corresponding to roof middle and floor middle
    ax.axvline(x=0.0, color='k', linestyle='--', linewidth=1)
    # ax.minorticks_on() ###to get subgrid
    # ax.grid(True)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlabel(r'Distance [m]', fontsize=21)
    if j==0:
        ax.set_ylabel(r'z [m]', fontsize=21)


if __name__ == "__main__":

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/GrandeMotteTignes/.')
    #########################################################################
    ####       BELOW CODE LINES TO PRODUCE FIGURES FOR ALL SIMUS         ####
    #########################################################################
    #### Transect and corresponding tunnel floor altitude (1% slope)
    Col_Names_Map = ['Node','X','Y','Z','Vx','Vy','Vz']
    Col_Names_Contour = ['Node','X','Y','Z']
    Transect = 'Transect4'
    Shape = 'Ovoide'
    Case = 'HalfChannelSyntheSym_NoSlid'
    tsp_list = [10,12]

    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###~~~~~~~~         GET DATA     ~~~~~~~~###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###################################################
    ########## START LOOP ON TIMESTEP #################
    ###################################################
    for k,tsp in enumerate(tsp_list):
        ###Get proper axe
        ax = axes[k]
        ###Set title of subplot
        ax.set_title('@{} days'.format(round(tsp*365/60)), fontsize=21, weight='bold')
        ax.set_aspect('equal', adjustable='box')  ###same horizontal and vertical scales
        ###Open the elmer output file
        filename_contour = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/ContourDomain_Node_X_Y_Z_{2}_60tspPerY_tsp{3}.dat'.format(Transect, Shape, Case, tsp)
        filename_map= 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/Node_X_Y_Z_U_V_W_{2}_60tspPerY_tsp{3}.dat'.format(Transect, Shape, Case, tsp)
        #####Load the Elmer output file
        df_contour = pd.read_csv(pathroot_mycode.joinpath(filename_contour), names=Col_Names_Contour,delim_whitespace=True)
        ### Close contour by copying first line and adding at the end of data frame
        df_contour_closed = pd.concat([df_contour, df_contour.iloc[[0]]], ignore_index=True)
        ###NOW LOAD DATA OVER WHOLE DOMAIN FOR MAP
        df_map = pd.read_csv(pathroot_mycode.joinpath(filename_map), names=Col_Names_Map,delim_whitespace=True)
        ###Create refined regular grid and interpolate field over grid
        x = np.arange(np.floor(np.min(df_map['X'])) - 1, np.floor(np.max(df_map['X'])) + 1, 0.1)
        y = np.arange(np.floor(np.min(df_map['Y'])) - 1, np.floor(np.max(df_map['Y'])) + 1, 0.1)
        X, Y = np.meshgrid(x, y)
        Vx_MeterPerYear = Interpolate_field(df_map, 'Vx', X, Y)
        Vx = 100/365*Vx_MeterPerYear ###convert in cm/day
        # shading
        clevs = np.arange(-1.6, 0.0, 0.1)  ## cbar for shading
        cmap = 'coolwarm'
        # colorbar
        levs_ticks = np.arange(-1.6, 0.2, 0.2)
        ###Fills up the map with colors for Vx
        CS1 = ax.contourf(X, Y, Vx, clevs, cmap=cmap, extend='both')
        ax.plot(df_contour_closed ['X'].values, df_contour_closed ['Y'].values, color='k', linewidth=2)
        ###Show colorbar (just once)
        if k ==1:
            cbar = plt.colorbar(CS1, ticks=levs_ticks, orientation='vertical', label=r'$V_\mathrm{x}$ [cm/day]')
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[df_contour_closed ['X'].values, df_contour_closed ['Y'].values])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in CS1.collections:
            c.set_clip_path(patch)



    plt.show()









