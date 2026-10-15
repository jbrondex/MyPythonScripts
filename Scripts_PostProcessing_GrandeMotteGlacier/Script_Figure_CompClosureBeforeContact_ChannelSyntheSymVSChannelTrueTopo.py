"""
@author: jbrondex

Description:
------------
This script is used to produce figure showing deformation of channel over time and more specifically to
compare contours obtained for the true topography versus the ones obtained when using the synthetic symmetric
counterpart before the contact between walls is established. We do that on a transect by transect basis.
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

###Fig 1 is subplots for right before and right after contact
fig, axes = plt.subplots(1, 3, figsize=(50, 50), sharey=True, constrained_layout=True)

if __name__ == "__main__":

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/GrandeMotteTignes/.')
    #########################################################################
    ####       BELOW CODE LINES TO PRODUCE FIGURES FOR ALL SIMUS         ####
    #########################################################################
    #### Transect and corresponding tunnel floor altitude (1% slope)
    Col_Names_Transient = ['Tsp','CallCount','BC','Node','X','Y','Z','Vx','Vy','SigmaI','Vn','Sigma_nn']
    Transect = 'Transect4'
    Floor_altitudes = [2803.5, 2803, 2802.5, 2802, 2801.5, 2801.12]
    BottomTunnel_altitudes = [2792.2, 2791.7, 2791.2, 2790.7, 2790.2, 2789.82]
    Bed_altitudes = [2784.75, 2784.75, 2779.53, 2774.56, 2773.12, 2774.82]
    #### Cases before/after incision
    Cases=['Channel_NoSlid']
    ###Times in months at which we want to plot the deformed tunnel (and corresponding names for title)
    Times = [1*5, 2*5, 9*5]##5 tsp each month for the case with 60 tsp per year
    Time_Names = ['@1m', '@4m', '@9m']
    ################################################
    ########## START LOOP ON CASES #################
    ################################################
    for Case in Cases: ###One Fig per case
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~~~~~~         GET TUNNEL CONTOUR     ~~~~~~~~###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###Get the bed and top DEM of the domain contour
        file_name_bed = './Maillages/DEMs/DEM_Bed_{}.dat'.format(Transect)
        # Load DEMs of domain
        Col_Names = ['x', 'y']
        DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
        ###Get ybed at x=25 for placing the cross
        df_ybed = DEM_bed.iloc[(DEM_bed['x'] - 25).abs().argsort()[:2]]
        ybedmiddle = (df_ybed.iloc[0]['y'] + df_ybed.iloc[1]['y']) / 2
        ###Open the elmer output file
        if Case == 'Channel_WithSlid':
            filename_output = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_HalfChannelSyntheSym_{0}_Ovo_Slid_60StepPerY_.dat'.format(Transect)
            filename_output_trueDEM = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_Channel_{0}_Ovo_Slid_60StepPerY_.dat'.format(Transect)
        elif Case == 'Channel_NoSlid':
            filename_output = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_HalfChannelSyntheSym_{0}_Ovo_NoSlid_60StepPerY_.dat'.format(Transect)
            filename_output_trueDEM = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_Channel_{0}_Ovo_NoSlid_60StepPerY_.dat'.format(Transect)
        #####Load the Elmer output file
        df0 = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Transient,delim_whitespace=True)
        df0_trueDEM = pd.read_csv(pathroot_mycode.joinpath(filename_output_trueDEM), names=Col_Names_Transient,delim_whitespace=True)
        ###################################################
        ########## START LOOP ON OUTPUT TIME ##############
        ###################################################
        for j,(Time,Time_Name) in enumerate(zip(Times,Time_Names)):
            ### Get proper subplot
            ###Get proper axe
            ax = axes[j]
            ###Set title of subplot
            ax.set_title('@{}m'.format(round(Time * 12 / 60)), fontsize=21, weight='bold')
            ax.set_aspect('equal', adjustable='box')  ###same horizontal and vertical scales
            ###Treatment of Data for the true DEM case
            df_trueDEM = df0_trueDEM[df0_trueDEM['Tsp'] == Time]
            ### get (closed) contour of tunnel for that time
            xt_init_trueDEM = df_trueDEM['X'].values
            xt_init_trueDEM = np.append(xt_init_trueDEM, xt_init_trueDEM[0])
            yt_init_trueDEM = df_trueDEM['Y'].values
            yt_init_trueDEM = np.append(yt_init_trueDEM, yt_init_trueDEM[0])
            ###Sort coordinates of contour clockwise
            xy_trueDEM = np.array([xt_init_trueDEM, yt_init_trueDEM])
            xy_trueDEM = np.transpose(xy_trueDEM)
            df_xy_sorted_trueDEM = pd.DataFrame({'X': xy_trueDEM[:, 0], 'Y': xy_trueDEM[:, 1]})
            ####Here we spot on paraview all cases that are closed at considered time to avoid plotting them
            Closed = False
            if (Case == 'Channel_WithSlid'):
                if (Transect in ['Transect3', 'Transect4', 'Transect5', 'Transect6']):
                    Closed = True
                elif Transect == 'Transect1':
                    if Time > 7*5:
                        Closed = True
                elif Transect == 'Transect2':
                    if Time > 3*5:
                        Closed = True
            elif (Case == 'Channel_NoSlid'):
                if Transect == 'Transect2':
                    if Time > 8*5:
                        Closed = True
                elif Transect == 'Transect3':
                    if Time > 3*5:
                        Closed = True
                elif Transect == 'Transect4' or Transect == 'Transect5':
                    if Time > 2*5:
                        Closed = True
                elif Transect == 'Transect6':
                    if Time > 2*5:
                        Closed = True
            #### NOW TREAT THE DATA FOR THE SYNTHETIC DOMAIN CASE #####
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
            ###Set the axe limits and ticks
            xmin = 17
            xmax = 33
            ###ymin is relative to bed altitude
            df_ybed = DEM_bed.iloc[(DEM_bed['x'] - 25).abs().argsort()[:2]]
            ymin = Bed_altitudes[int(Transect[-1])-1] - 1
            ymax = 2807
            ax.set_xticks([19, 25, 31])
            ax.set_yticks([2775, 2780, 2785, 2790, 2795, 2800, 2805])
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.tick_params(labelsize=16)  # fontsize of the tick labels
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.set_aspect('equal')
            ax.set_xlabel(r'Distance [m]', fontsize=21)
            if j == 0:
                ax.set_ylabel(r'z [m]', fontsize=21)
            # ax.grid(True)
            ###Plot the stuff###
            ###print contour of domain
            ax.plot(DEM_bed['x'].values, DEM_bed['y'].values, color='k', linewidth=2)
            ###Brown color for bed
            ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float),2750, color='tan')
            ###Icy color for ice
            ax.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float),2830, color='lightblue')
            ###Plot contour of tunnel and fill in in white if tunnel is not closed
            if not (Closed):
                ax.add_patch(Patch.Polygon(np.transpose(np.array([df_xy_sorted_trueDEM['X'].values, df_xy_sorted_trueDEM['Y'].values])),closed=True, fill=True, facecolor='w', edgecolor='k'))
                ###Add contour of synthetic symmetric case
                if not WallTouching:
                    ax.plot(xfull+25, yfull, color='r', linestyle='--', linewidth=1.8)
                else:
                    ax.plot(xupperfull+25, yupperfull, color='r', linestyle='--', linewidth=1.8)
                    ax.plot(xlowerfull+25, ylowerfull, color='r', linestyle='--', linewidth=1.8)
            else:  ##otherwise plot a cross
                ax.scatter(0.5 * (xmin + xmax), 0.5 * (ybedmiddle + ymax), marker='x', s=80, c='k')



    plt.show()


        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~         SAVE THE PLOTS (ONE PER SHAPE/WIDTH)          ~~###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # ###Save figure with full transect profile and initial tunnel shape
        # name_output_fig4 = 'DeformationAt1_4_9months_AllTransects_{0}_HalfChannelSyntheSym'.format(Case)
        # name_path_fig4 = '/home/brondexj/BETTIK/GrandeMotteTignes/RunTransient_Closure/Postprocessing/Figures/.'
        # path_output_fig4 = Path(name_path_fig4)
        # fig4.savefig(path_output_fig4.joinpath(name_output_fig4))






