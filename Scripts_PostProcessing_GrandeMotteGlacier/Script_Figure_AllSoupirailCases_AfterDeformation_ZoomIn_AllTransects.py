"""
@author: jbrondex

Description:
------------
This script is used to produce figure showing deformation of soupirail over time. Beware that the tunnel cases (init, incised straight walls,
incised incurved walls) are treated in a separated script. A figure is made of six subplots (one per transect). Each subplots is divided in three:
@6months, @9months, @12months. We produce one figure for each case (with sliding, without sliding) but for ovoid shape only.
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
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#######We define a function to interpolate bottom corners of soupirail on bed DEM
def Interpolate_on_bed(DEM_bed, xcorner):
    xbed = DEM_bed['x']
   # xbed = xbed[::-1]  ###the scipy interpolation function requires the x to be sorted in increasing order
    ybed = DEM_bed['y']
    #ybed = ybed[::-1]  ###so that xbed and ybed are consistent
    # interpolation = interpolate.CubicSpline(xbed, ybed)
    # ycorner = interpolation(xcorner)
    ycorner = np.interp(xcorner,xbed,ybed)
    return ycorner

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
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=24)



if __name__ == "__main__":

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/GrandeMotteTignes/.')
    #########################################################################
    ####       BELOW CODE LINES TO PRODUCE FIGURES FOR ALL SIMUS         ####
    #########################################################################
    #### Transect and corresponding tunnel floor altitude (1% slope)
    Col_Names_Transient = ['Tsp','CallCount','BC','Node','X','Y','Z','Vx','Vy','SigmaI','Vn','Sigma_nn']
    Transect_Names = ['Transect1', 'Transect2', 'Transect3', 'Transect4', 'Transect5', 'Transect6']
    Floor_altitudes = [2803.5, 2803, 2802.5, 2802, 2801.5, 2801.12]
    BottomTunnel_altitudes = [2792.2, 2791.7, 2791.2, 2790.7, 2790.2, 2789.82]
    #### Cases before/after incision
    Cases=['Soupirail_WithSlid', 'Soupirail_NoSlid']
    #### Possible shapes of the tunnel
    Shapes = ['Ovoide']
    #### Width of tunnel (for Rectangle and Half Circle case only) and corresponding name (in cm) for output file name
    Widths = [1, 1.5, 2]
    Width_Names = ['W100cm','W150cm', 'W200cm']
    ###Times in months at which we want to plot the deformed tunnel (and corresponding names for title)
    Times = [6, 9, 12]
    Time_Names = ['@6m', '@9m', '@1a']
    ######################################################
    ########## START LOOP ON SHAPES ######################
    ######################################################
    for Shape in Shapes: ###One Figure per shape
        ######################################################
        ########## START LOOP ON WIDTHS ######################
        ######################################################
        for l, (Width, Width_Name) in enumerate(zip(Widths, Width_Names)):
            if Shape == 'Ovoide' and l>0:
                break###Only one width for ovoide shape
            ################################################
            ########## START LOOP ON CASES #################
            ################################################
            for Case in Cases: ###One Fig per case
                print('Dealing with Figure for ', Shape, ' with width ', Width, ' and ', Case)
                ###Parameters of figure
                fig4 = plt.figure(figsize=(16, 8))
                outer = gridspec.GridSpec(2, 3, wspace=0.2, hspace=0.05)  ##9 Main Subplots, each divided in 2
                fig4.text(0.5, 0.06, r'Distance [m]', ha='center', fontsize=22)
                fig4.text(0.05, 0.5, r'z [m]', va='center', rotation='vertical', fontsize=22)
                ################################################
                ########## START LOOP ON TRANSECTS##############
                ################################################
                for k, Transect in enumerate(Transect_Names):
                    ### Get proper main subplot (one per transect
                    axes = np.empty(shape=(1, 3), dtype=object)  ##This will be for subsubplots
                    inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[k], hspace=0, wspace=0)
                    ###Get the bed and top DEM of the domain contour
                    file_name_bed = './Maillages/DEMs/DEM_Bed_{}.dat'.format(Transect)
                    # Load DEMs of domain
                    Col_Names = ['x', 'y']
                    DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names,delim_whitespace=True)
                    ###Get ybed at x=25 for placing the cross
                    df_ybed = DEM_bed.iloc[(DEM_bed['x'] - 25).abs().argsort()[:2]]
                    ybedmiddle = Interpolate_on_bed(DEM_bed, 25)

                    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                    ###~~~~~~~~         GET TUNNEL CONTOUR     ~~~~~~~~###
                    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                    ###Open the elmer output file
                    if Shape == 'Ovoide':
                        if Case == 'Soupirail_WithSlid':
                            filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_Soupirail_{0}_Ovo_Slid_.dat'.format(Transect, Shape)
                        elif Case == 'Soupirail_NoSlid':
                            filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_Soupirail_{0}_Ovo_NoSlid_.dat'.format(Transect, Shape)
                    elif Shape == 'Rectangle':
                        if Case == 'Soupirail_WithSlid':
                            filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_Soupirail_{0}_Rect_{2}_Slid_.dat'.format(Transect, Shape, Width_Name)
                        elif Case == 'Soupirail_NoSlid':
                            filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_Soupirail_{0}_Rect_{2}_NoSlid_.dat'.format(Transect, Shape, Width_Name)
                    elif Shape == 'Circle':
                        if Case == 'Soupirail_WithSlid':
                            filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_Soupirail_{0}_Circ_{2}_Slid_.dat'.format(Transect, Shape, Width_Name)
                        elif Case == 'Soupirail_NoSlid':
                            filename_output = 'RunTransient_Closure/RunTransient_{0}_{1}/ScalarOutput/TunnelOutput_Soupirail_{0}_Circ_{2}_NoSlid_.dat'.format(Transect, Shape, Width_Name)
                    #####Load the Elmer output file
                    df0 = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Transient,delim_whitespace=True)
                    ###################################################
                    ########## START LOOP ON OUTPUT TIME ##############
                    ###################################################
                    for j,(Time,Time_Name) in enumerate(zip(Times,Time_Names)):
                        ###Get Elmer output corresponding to considered time
                        df = df0[df0['Tsp']==Time]
                        ### get (closed) contour of tunnel for that time
                        xt_init = df['X'].values
                        yt_init = df['Y'].values
                        ### for transect 4 remove righter-most points
                        if Transect == 'Transect4':
                            xt_init=xt_init[2::]
                            yt_init=yt_init[2::]
                            ###Bring back rightmost point on bed
                            yt_right_on_bed = Interpolate_on_bed(DEM_bed, xt_init[0])
                            yt_init[0]=yt_right_on_bed
                        elif Transect == 'Transect5':##Same for Transect 5 but for only one point
                            xt_init=xt_init[1::]
                            yt_init=yt_init[1::]
                            ###Bring back rightmost point on bed
                            yt_right_on_bed = Interpolate_on_bed(DEM_bed, xt_init[0])
                            yt_init[0]=yt_right_on_bed
                        elif Transect == 'Transect6':##For T6 just reproject first and last point on bed
                            yt_right_on_bed = Interpolate_on_bed(DEM_bed, xt_init[0])
                            yt_init[0]=yt_right_on_bed
                            yt_left_on_bed = Interpolate_on_bed(DEM_bed, xt_init[-1])
                            yt_init[-1]=yt_left_on_bed

                        # yt_init[-1] = Interpolate_on_bed(DEM_bed,xt_init[-1])
                        ###Add middle point of bed
                        # xt_init = np.append(xt_init, 25)
                        # yt_init = np.append(yt_init, ybedmiddle)
                        ###close de contour by re-adding first point
                        xt_init = np.append(xt_init, xt_init[0])
                        yt_init = np.append(yt_init, yt_init[0])
                        ##Sort coordinates of contour clockwise
                        xy = np.array([xt_init, yt_init])
                        xy = np.transpose(xy)
                        df_xy_sorted = pd.DataFrame({'X': xy[:, 0], 'Y': xy[:, 1]})
                        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                        ###~~~         GET PROPER SUB-SUBPLOT AND PLOT          ~~###
                        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                        ###Get proper sub-subplots
                        ax = plt.Subplot(fig4, inner[j])
                        ###Name of subsubplots
                        ax.set_title('{}'.format(Time_Name), fontsize=16)
                        ###Only show axis that are needed
                        if k <3:
                            ax.xaxis.set_visible(False) ###xaxis only for bottom row
                        if j > 0:
                            ax.yaxis.set_visible(False) ###no yaxis between subsubplots
                        ###Add subtitles to subsubplots corresponding to years
                        # Add ghost axes and titles on columns
                        ax_transect = fig4.add_subplot(inner[:])
                        ax_transect.axis('off')
                        ax_transect.set_title('Transect {}'.format(k+1), fontsize=20, weight='bold')
                        ###Set the axe limits and ticks
                        xmin=23.7
                        xmax=26.3
                        ###ymin is relative to bed altitude
                        df_ybed = DEM_bed.iloc[(DEM_bed['x'] - 25).abs().argsort()[:2]]
                        ymin=ybedmiddle-1.5
                        ymax=ybedmiddle+2.5
                        ax.set_xticks([24,25,26])
                        ax.set_yticks(list(range(2770,2820,2)))
                        ax.set_xlim(xmin, xmax)
                        ax.set_ylim(ymin, ymax)
                        ax.tick_params(labelsize=16)  # fontsize of the tick labels
                        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                        ax.set_aspect('equal')
                        # ax.grid(True)
                        axes[0, j] = ax
                        ###Plot the stuff###
                        ###Brown color for bed
                        ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2750, color='tan')
                        ###Icy color for ice
                        ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2830, color='lightblue')
                        ###Plot contour of soupirail and fill in in white
                        ax.add_patch(Patch.Polygon(np.transpose(np.array([df_xy_sorted['X'].values,df_xy_sorted['Y'].values])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                        ###print contour of domain
                        ax.plot(DEM_bed['x'].values, DEM_bed['y'].values, color='k', linewidth=3)
                        fig4.add_subplot(ax)

                plt.show()


                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###~~~         SAVE THE PLOTS (ONE PER SHAPE/WIDTH)          ~~###
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###Save figure with full transect profile and initial tunnel shape
                name_output_fig4 = 'DeformationAt6_9_12months_AllTransects_{0}_Ovoide'.format(Case,Shape)
                name_path_fig4 = '/home/brondexj/BETTIK/GrandeMotteTignes/RunTransient_Closure/Postprocessing/Figures/.'
                path_output_fig4 = Path(name_path_fig4)
                fig4.savefig(path_output_fig4.joinpath(name_output_fig4))

                plt.show()




