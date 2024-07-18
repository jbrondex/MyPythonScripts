"""
@author: jbrondex

Description:
------------
This scripts is used to produce figure showing a zoom-in on , for the three shapes, the tunnel init, the tunnel incised, the channel. This is done for Transect1 only
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

###Parameters of figure
fig2, axes = plt.subplots(4, 3, figsize=(12, 17), sharey= True)
fig2.subplots_adjust(hspace=0.11,wspace=0.0)
for i in range(0,4):
    for j in range(0,3):
        ax = axes[i,j]
        ax.set_xlim([15, 35])
        ax.set_aspect('equal')
        ax.locator_params(axis = 'y', nbins = 4)
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        # ax.minorticks_on() ###to get subgrid
        ax.grid(True)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        if j==0:
            ax.set_ylabel(r'z [m]', fontsize=21)
        if i==3:
            ax.set_xlabel(r'Distance [m]', fontsize=21)

if __name__ == "__main__":

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/GrandeMotteTignes/.')
    #########################################################################
    ####       BELOW CODE LINES TO PRODUCE FIGURES FOR ALL SIMUS         ####
    #########################################################################
    #### Transect and corresponding tunnel floor altitude (1% slope)
    Col_Names_Transient = ['Tsp','CallCount','BC','Node','X','Y','Z','Vx','Vy','SigmaI','Vn','Sigma_nn']
    Transect = 'Transect1'
    Floor_altitudes = [2803.5]
    #### Cases before/after incision
    Cases=['TunnelInit', 'TunnelIncised_StraightWalls', 'TunnelIncised_IncurvedWalls', 'Channel']
    #### Possible shapes of the tunnel
    Shapes = ['Ovoide', 'Circle', 'Rectangle']
    #### Width of tunnel (for Rectangle and Half Circle case only) and corresponding name (in cm) for output file name
    Widths = [1, 1.5, 2]
    Width_Names = ['W100cm','W150cm', 'W200cm']
    Width_Name = 'W200cm'
    ######################################################
    ########## START LOOP ON SHAPES ######################
    ######################################################
    for j,Shape in enumerate(Shapes): ###One column per shape
        ################################################
        ########## START LOOP ON CASES #################
        ################################################
        for i, Case in enumerate(Cases): ###One line per case
            print('Transect n°', Transect, '   Shape:', Shape, '   Case:', Case)
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
                if Case == 'Channel':
                    filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Ovo_NoSlid_.dat'.format(Transect, Shape, Case)
                else:
                    filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Ovo_.dat'.format(Transect, Shape, Case)
            elif Shape == 'Rectangle':
                if Case == 'Channel':
                    filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_{3}_NoSlid_.dat'.format(Transect, Shape, Case, Width_Name)
                    # filename_output_w100 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_W100cm_NoSlid_.dat'.format(Transect, Shape, Case)
                    # filename_output_w150 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_W150cm_NoSlid_.dat'.format(Transect, Shape, Case)
                    # filename_output_w200 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_W200cm_NoSlid_.dat'.format(Transect, Shape, Case)
                else:
                    filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_{3}_.dat'.format(Transect, Shape, Case, Width_Name)
                    # filename_output_w100 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_W100cm_.dat'.format(Transect, Shape, Case)
                    # filename_output_w150 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_W150cm_.dat'.format(Transect, Shape, Case)
                    # filename_output_w200 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Rect_W200cm_.dat'.format(Transect, Shape, Case)
            elif Shape == 'Circle':
                if Case == 'Channel':
                    filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_{3}_NoSlid_.dat'.format(Transect, Shape, Case, Width_Name)
                    # filename_output_w100 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_W100cm_NoSlid_.dat'.format(Transect, Shape, Case)
                    # filename_output_w150 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_W150cm_NoSlid_.dat'.format(Transect, Shape, Case)
                    # filename_output_w200 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_W200cm_NoSlid_.dat'.format(Transect, Shape, Case)
                else:
                    filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_{3}_.dat'.format(Transect, Shape, Case, Width_Name)
                    # filename_output_w100 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_W100cm_.dat'.format(Transect, Shape, Case)
                    # filename_output_w150 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_W150cm_.dat'.format(Transect, Shape, Case)
                    # filename_output_w200 = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_{2}_{0}_Circ_W200cm_.dat'.format(Transect, Shape, Case)
            df0 = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Transient,delim_whitespace=True)
            # if not (Shape=='Ovoide'):
            #     df0_w100 = pd.read_csv(pathroot_mycode.joinpath(filename_output_w100), names=Col_Names_Transient,delim_whitespace=True)
            #     df0_w150 = pd.read_csv(pathroot_mycode.joinpath(filename_output_w150), names=Col_Names_Transient,delim_whitespace=True)
            #     df0_w200 = pd.read_csv(pathroot_mycode.joinpath(filename_output_w200), names=Col_Names_Transient,delim_whitespace=True)
            #### Keep only first tsp
            df = df0[df0['Tsp']==1]
            # if not (Shape=='Ovoide'):
            #     df_w100 = df0_w100[df0_w100['Tsp']==1]
            #     df_w150  = df0_w150 [df0_w150['Tsp']==1]
            #     df_w200  = df0_w200 [df0_w200 ['Tsp']==1]
            ### get (closed) contour of initial tunnel
            xt_init = df['X'].values
            xt_init = np.append(xt_init, xt_init[0])
            yt_init = df['Y'].values
            yt_init = np.append(yt_init, yt_init[0])
            ###Sort coordinates of contour clockwise
            xy = np.array([xt_init, yt_init])
            xy = np.transpose(xy)
            xy_sorted = sort_coordinates(xy)
            # if not (Shape=='Ovoide'):
            #     ###case w100
            #     xt_init_w100 = df_w100['X'].values
            #     xt_init_w100 = np.append(xt_init_w100, xt_init_w100[0])
            #     yt_init_w100 = df_w100['Y'].values
            #     yt_init_w100 = np.append(yt_init_w100, yt_init_w100[0])
            #     ###Sort coordinates of contour clockwise
            #     xy_w100 = np.array([xt_init_w100, yt_init_w100])
            #     xy_w100 = np.transpose(xy_w100)
            #     xy_w100_sorted = sort_coordinates(xy_w100)
            #     ###case w150
            #     xt_init_w150 = df_w150['X'].values
            #     xt_init_w150 = np.append(xt_init_w150, xt_init_w150[0])
            #     yt_init_w150 = df_w150['Y'].values
            #     yt_init_w150 = np.append(yt_init_w150, yt_init_w150[0])
            #     ###Sort coordinates of contour clockwise
            #     xy_w150 = np.array([xt_init_w150, yt_init_w150])
            #     xy_w150 = np.transpose(xy_w150)
            #     xy_w150_sorted = sort_coordinates(xy_w150)
            #     ###case w200
            #     xt_init_w200 = df_w200['X'].values
            #     xt_init_w200 = np.append(xt_init_w200, xt_init_w200[0])
            #     yt_init_w200 = df_w200['Y'].values
            #     yt_init_w200 = np.append(yt_init_w200, yt_init_w200[0])
            #     ###Sort coordinates of contour clockwise
            #     xy_w200 = np.array([xt_init_w200, yt_init_w200])
            #     xy_w200 = np.transpose(xy_w200)
            #     xy_w200_sorted = sort_coordinates(xy_w200)

            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###~~~         PREPARE AND MAKE THE PLOTS          ~~###
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###Get proper axe
            ax = axes[i, j]
            ###Set title of subplot
            # ax.set_title('Transect {}'.format(k + 1), fontsize=21, weight='bold')
            # ax.set_aspect('equal', adjustable='box') ###same horizontal and vertical scales
            ###print contour of domain
            ax.plot(xc, yc, color='k', linewidth=2)
            ###Brown color for bed
            ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2750, color='tan')
            ###Icy color for ice
            ax.fill_between(np.array(DEM_top['x'].values, dtype=float), np.array(DEM_top['y'].values, dtype=float),2799, color='lightblue')
            ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2801, color='lightblue')
            ###Plot contour of tunnel and fill in in white
            ax.add_patch(Patch.Polygon(xy_sorted, closed=True, fill=True, facecolor= 'w', edgecolor='k' ))
            ###For the circle and rectangle shape plot the width 100cm and 150cm cases
            # if not (Shape=='Ovoide'):
            #     ax.plot(xy_w100_sorted[0],xy_w100_sorted[1],color='k', linewidth=1, linestyle=':')
            #     ax.plot(xy_w150_sorted[0],xy_w150_sorted[1],color='k', linewidth=1, linestyle='-')
            ###Fix ylim depending on case
            ymin=np.min(yt_init)-1
            ymax=np.max(yt_init)+1
            ax.set_ylim([ymin, ymax])

    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###~~~         SAVE THE PLOTS (ONE PER SHAPE/WIDTH)          ~~###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    name_output_fig2 = 'AllCases_AllShapes_@0a_{}_'.format(Transect)
    name_path_fig2 = '/home/brondexj/BETTIK/GrandeMotteTignes/RunTransient_Closure/Postprocessing/Figures/.'
    path_output_fig2 = Path(name_path_fig2)
    fig2.savefig(path_output_fig2.joinpath(name_output_fig2))

    plt.show()