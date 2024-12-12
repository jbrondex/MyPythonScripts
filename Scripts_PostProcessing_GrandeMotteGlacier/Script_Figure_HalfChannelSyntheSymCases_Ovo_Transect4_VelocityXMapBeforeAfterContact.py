"""
@author: jbrondex

Description:
------------
This script is used to produce figure showing deformation of channel over time. Beware that the tunnel cases (init, incised straight walls,
incised incurved walls) are treated in a separated script. A figure is made of six subplots (one per transect). Each subplots is divided in three:
@ different times. We produce one figure for each shape and for each case (with sliding, without sliding)
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
    Bed_altitudes = [2784.75, 2784.75, 2779.53, 2774.56, 2773.12, 2774.82]
    #### Cases before/after incision
    Cases=['Channel_WithSlidArg']#, 'Channel_NoSlid']
    ###Times in months at which we want to plot the deformed tunnel (and corresponding names for title)
    Times = [1*5, 4*5, 9*5]##5 tsp each month for the case with 60 tsp per year
    Time_Names = ['@1m', '@4m', '@9m']
    ################################################
    ########## START LOOP ON CASES #################
    ################################################
    for Case in Cases: ###One Fig per case
        print('Dealing with Figure for ', Case)
        ###Parameters of figure
        fig4 = plt.figure(figsize=(16, 8))
        outer = gridspec.GridSpec(2, 3, wspace=0.15, hspace=0.05)  ##9 Main Subplots, each divided in 2
        fig4.text(0.5, 0.06, r'Distance [m]', ha='center', fontsize=22)
        fig4.text(0.05, 0.5, r'z [m]', va='center', rotation='vertical', fontsize=22)
        ################################################
        ########## START LOOP ON TRANSECTS##############
        ################################################
        for k, (Transect, Bed_altitude) in enumerate(zip(Transect_Names,Bed_altitudes)):
            ### Get proper main subplot (one per transect
            axes = np.empty(shape=(1, 3), dtype=object)  ##This will be for subsubplots
            inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[k], hspace=0, wspace=0)
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###~~~~~~~~         GET TUNNEL CONTOUR     ~~~~~~~~###
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###Open the elmer output file
            if Case == 'Channel_WithSlidArg':
                filename_output = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_HalfChannelSyntheSym_{0}_Ovo_SlidArg_60StepPerY_.dat'.format(Transect)
            elif Case == 'Channel_NoSlid':
                filename_output = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_HalfChannelSyntheSym_{0}_Ovo_NoSlid_60StepPerY_.dat'.format(Transect)
            #####Load the Elmer output file
            df0 = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Transient,delim_whitespace=True)
            ###################################################
            ########## START LOOP ON OUTPUT TIME ##############
            ###################################################
            for j,(Time,Time_Name) in enumerate(zip(Times,Time_Names)):
                ###Get Elmer output corresponding to considered time
                dfright = df0[df0['Tsp']==Time]
                dfright.reset_index(inplace=True)
                ### check if contour contains negative x (meaning that walls are touching each other
                WallTouching = (dfright['X'].values < 0).any()
                if WallTouching:### If wall touching, then channel is divided into two closed contour that we want to distinguish
                    ###find index of first node (going from top to bottom) beyond symmetry axis
                    index_min=min(dfright.index[dfright['X'] < 0])
                    dfupperright=dfright.loc[0:index_min-1]
                    dfupperleft = dfupperright[1:] ##we remove the first node as it is on symmetry axis
                    ###Reverse dfleft to have nodes from bottom to top (clockwise contour)
                    dfupperleft = dfupperleft.iloc[::-1]
                    ###Get node coordinates of contour
                    xupperright = dfupperright['X'].values
                    xupperleft = -dfupperleft['X'].values
                    xupperfull = np.append(xupperright, xupperleft)
                    yupperright = dfupperright['Y'].values
                    yupperleft = dfupperleft['Y'].values
                    yupperfull = np.append(yupperright, yupperleft)
                    ###NOW SAME THING FOR BOTTOM PART OF CONTOUR###
                    ###find index of last node beyond symmetry axis
                    index_max=max(dfright.index[dfright['X'] < 0])
                    dflowerright=dfright.loc[index_max+1::]
                    dflowerleft = dflowerright[1:] ##we remove the first node as it is on symmetry axis
                    ###Reverse dfleft to have nodes from bottom to top (clockwise contour)
                    dflowerleft = dflowerleft.iloc[::-1]
                    ###Get node coordinates of contour
                    xlowerright = dflowerright['X'].values
                    xlowerleft = -dflowerleft['X'].values
                    xlowerfull = np.append(xlowerright, xlowerleft)
                    ylowerright = dflowerright['Y'].values
                    ylowerleft = dflowerleft['Y'].values
                    ylowerfull = np.append(ylowerright, ylowerleft)
                else:##No wall touching so one single closed contour built from main data frame
                    ###Build the symmetric wall
                    dfleft = dfright[1:] ##we remove the first node as it is on symmetry axis
                    ###Reverse dfleft to have nodes from bottom to top (clockwise contour)
                    dfleft = dfleft.iloc[::-1]
                    ###Get node coordinates of contour
                    xright=dfright['X'].values
                    xleft=-dfleft['X'].values
                    xfull=np.append(xright,xleft)
                    yright=dfright['Y'].values
                    yleft=dfleft['Y'].values
                    yfull=np.append(yright,yleft)



                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###~~~         GET PROPER SUB-SUBPLOT AND PLOT          ~~###
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
                ###Get proper sub-subplots
                ax = plt.Subplot(fig4, inner[j], sharey=axes[0, 0])
                ###Name of subsubplots
                ax.set_title('{}'.format(Time_Name), fontsize=16)
                ###Only show axis that are needed
                if not (k==0 or k==3): ##yaxis only for first column
                    ax.yaxis.set_visible(False)
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
                xmin=-10
                xmax=10
                ###ymin is relative to bed altitude
                ymin=2773
                ymax=2806
                ax.set_xticks([-6,0,6])
                ax.set_yticks([2775,2780,2785,2790,2795,2800,2805])
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.tick_params(labelsize=16)  # fontsize of the tick labels
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                ax.set_aspect('equal')
                # ax.grid(True)
                axes[0, j] = ax
                ###Plot the stuff###
                ###print contour of domain
                ax.axhline(y=Bed_altitude, color='k', linestyle='-', linewidth=2)
                ###Brown color for bed
                ax.fill_between( [xmin,xmax], [Bed_altitude,Bed_altitude], 2750, color='tan')
                ###Icy color for ice
                ax.fill_between([xmin,xmax], [Bed_altitude,Bed_altitude], 2830, color='lightblue')
                ###Plot contour of tunnel and fill in in white if tunnel is not closed
                if WallTouching:
                    ax.add_patch(Patch.Polygon(np.transpose(np.array([xupperfull,yupperfull])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                    ax.add_patch(Patch.Polygon(np.transpose(np.array([xlowerfull,ylowerfull])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                else:
                    ax.add_patch(Patch.Polygon(np.transpose(np.array([xfull,yfull])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                fig4.add_subplot(ax)

        plt.show()


        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~         SAVE THE PLOTS (ONE PER SHAPE/WIDTH)          ~~###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###Save figure with full transect profile and initial tunnel shape
        name_output_fig4 = 'DeformationAt1_4_9months_AllTransects_{0}_HalfChannelSyntheSym'.format(Case)
        name_path_fig4 = '/home/brondexj/BETTIK/GrandeMotteTignes/RunTransient_Closure/Postprocessing/Figures/.'
        path_output_fig4 = Path(name_path_fig4)
        fig4.savefig(path_output_fig4.joinpath(name_output_fig4))






