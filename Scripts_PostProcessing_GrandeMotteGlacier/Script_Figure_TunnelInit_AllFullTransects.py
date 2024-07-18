"""
@author: jbrondex

Description:
------------
This scripts is used to produce figure showing the six whole transects as subplots with the tunnel in its initial shape for the 3 shapes and the 3 widths (in the
case of rectangle/half-circle shapes)
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
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y)
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

##Function below is used to sort point of contour clockwise
def sort_coordinates(list_of_xy_coords):
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
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

###Fig 1 is subplots for the steady max tensile principle stress around the initial tunnel for the 6 transects
###Parameters if figure not saved
fig1, axes = plt.subplots(2, 3, figsize=(15, 15), sharey=True)
fig1.subplots_adjust(hspace=0.09,wspace=0.1)
for i in range(0,2):
    for j in range(0,3):
        ax = axes[i,j]
        ax.set_xlim([0.0, 50.0])
        ax.set_ylim([2770, 2835])
        ax.locator_params(axis = 'y', nbins = 10)
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        # ax.minorticks_on() ###to get subgrid
        ax.grid(True)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        if i==0: ##Top Line Horizontal profile
            if j==0:
                ax.set_ylabel(r'z [m]', fontsize=21)
        elif i==1: ##Bottom Line Vertical profile
            # ax.set_ylim([0.0, 50.0])
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
    Col_Names_Steady = ['Tsp','BC','Node','X','Y','Z','Vx','Vy','SigmaI','Vn','Sigma_nn']
    Transect_Names = ['Transect1', 'Transect2', 'Transect3', 'Transect4', 'Transect5', 'Transect6']
    Floor_altitudes = [2803.5, 2803, 2802.5, 2802, 2801.5, 2801.12]
    #### Possible shapes of the tunnel
    Shapes = ['Rectangle', 'Circle', 'Ovoide']
    #### Width of tunnel (for Rectangle and Half Circle case only) and corresponding name (in cm) for output file name
    Widths = [1, 1.5, 2]
    Width_Names = ['W100cm','W150cm', 'W200cm']
    ######################################################
    ########## START LOOP ON SHAPES ######################
    ######################################################
    for Shape in Shapes:
        ###################################################
        ########## START LOOP ON WIDTHS ###################
        ###################################################
        for k, (Width, Width_Name) in enumerate(zip(Widths, Width_Names)):
            if Shape == 'Ovoide' and k>0:
                break ###Only one width for ovoide shape
            ###################################################
            ########## START LOOP ON TRANSECTS ################
            ###################################################
            for k, (Transect, Floor_altitude) in enumerate(zip(Transect_Names, Floor_altitudes)):
                print('Transect n°', Transect, '   Shape:', Shape, '   Width:', Width_Name)
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###~~~         CONSTRUCT DOMAIN CONTOUR          ~~###
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###Get the bed and top DEM of the domain contour
                file_name_top = './Maillages/DEMs/DEM_Top_{}.dat'.format(Transect)
                file_name_bed = './Maillages/DEMs/DEM_Bed_{}.dat'.format(Transect)
                # Load DEMs of domain
                Col_Names = ['x', 'y']
                ###DEM top
                DEM_top = pd.read_csv(pathroot_mycode.joinpath(file_name_top), names=Col_Names, delim_whitespace=True)
                ###Add a last point at exaclty x=50
                interpolation = interpolate.CubicSpline(DEM_top['x'], DEM_top['y'])
                y50 = interpolation(50.0)
                DEM_top.loc[-1] = [50.0, y50]  # adding a row
                DEM_top = DEM_top.reset_index(drop=True)
                Npttop = len(DEM_top.index)
                ###DEM bed
                DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
                ###Add a last point at exaclty x=50
                interpolation = interpolate.CubicSpline(DEM_bed['x'], DEM_bed['y'])
                y50 = interpolation(50.0)
                DEM_bed.loc[-1] = [50.0, y50]  # adding a row
                ###Reverse DEM_bed dataframe so that nodes are consecutive
                DEM_bed = DEM_bed[::-1]
                DEM_bed = DEM_bed.reset_index(drop=True)
                Nptbed = len(DEM_bed.index)
                ###Merge both DEM to get full contour of domain
                cont = [DEM_top, DEM_bed]
                MainContour = pd.concat(cont)
                MainContour = MainContour.reset_index(drop=True)
                ###close the contour of domain
                xc = MainContour['x'].values
                xc = np.append(xc, xc[0])
                yc = MainContour['y'].values
                yc = np.append(yc, yc[0])
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###~~~         CONSTRUCT TUNNEL CONTOUR          ~~###
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###Open the elmer output file
                if Shape == 'Ovoide':
                    filename_output = 'RunSteady_Tunnel/RunSteady_{0}_{1}/ScalarOutput/TunnelOutput_RunSteady_Tunnel_{0}_{1}_TunnelFreeSurf_.dat'.format(Transect, Shape)
                else:
                    filename_output = 'RunSteady_Tunnel/RunSteady_{0}_{1}_{2}/ScalarOutput/TunnelOutput_RunSteady_Tunnel_{0}_{1}_{2}_TunnelFreeSurf_.dat'.format(Transect, Shape, Width_Name)
                df = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Steady,delim_whitespace=True)
                ### get (closed) contour of initial tunnel
                xt_init = df['X'].values
                xt_init = np.append(xt_init, xt_init[0])
                yt_init = df['Y'].values
                yt_init = np.append(yt_init, yt_init[0])
                ###Sort coordinates of contour clockwise
                xy = np.array([xt_init, yt_init])
                xy = np.transpose(xy)
                xy_sorted = sort_coordinates(xy)
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###~~~         PREPARE AND MAKE THE PLOTS          ~~###
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###Get proper axe
                if k <= 2:
                    i = 0
                    j = k
                else:
                    i = 1
                    j = k - 3
                ax = axes[i, j]
                ###Set title of subplot
                ax.set_title('Transect {}'.format(k + 1), fontsize=21, weight='bold')
                ax.set_aspect('equal', adjustable='box') ###same horizontal and vertical scales
                ###print contour of domain
                ax.plot(xc, yc, color='k', linewidth=2)
                ###Brown color for bed
                ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2750, color='tan')
                ###Icy color for ice
                ax.fill_between(np.array(DEM_top['x'].values, dtype=float), np.array(DEM_top['y'].values, dtype=float),2799, color='lightblue')
                ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2801, color='lightblue')
                ###Plot contour of tunnel and fill in in white
                ax.add_patch(Patch.Polygon(xy_sorted, closed=True, fill=True, facecolor= 'w', edgecolor='k' ))
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###~~~         SAVE THE PLOTS (ONE PER SHAPE/WIDTH)          ~~###
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###Save figure with full transect profile and initial tunnel shape
            if Shape=='Ovoide':
                name_output_fig1 = 'Full6Transects_TunnelInit_{}'.format(Shape)
            else:
                name_output_fig1 = 'Full6Transects_TunnelInit_{}_{}'.format(Shape,  Width_Name)
            name_path_fig1 = '/home/brondexj/BETTIK/GrandeMotteTignes/RunTransient_Closure/Postprocessing/Figures/.'
            path_output_fig1 = Path(name_path_fig1)
            fig1.savefig(path_output_fig1.joinpath(name_output_fig1))

    plt.show()