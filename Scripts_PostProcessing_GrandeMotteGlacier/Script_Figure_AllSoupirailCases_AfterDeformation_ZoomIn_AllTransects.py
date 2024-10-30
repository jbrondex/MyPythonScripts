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
    ###Times in months at which we want to plot the deformed tunnel (and corresponding names for title)
    Times = [6*5, 9*5, 12*5]##5 tsp each month for the case with 60 tsp per year
    Time_Names = ['@6m', '@9m', '@1a']
    ###Parameters of figure 1 (both slid and no slid on same fig)
    fig1 = plt.figure(figsize=(16, 8))
    outer1 = gridspec.GridSpec(2, 3, wspace=0.2, hspace=0.05)  ##9 Main Subplots, each divided in 2
    fig1.text(0.5, 0.06, r'Distance [m]', ha='center', fontsize=22)
    fig1.text(0.05, 0.5, r'z [m]', va='center', rotation='vertical', fontsize=22)
    ###Parameters of figure 2 (only slid case)
    fig2 = plt.figure(figsize=(16, 8))
    outer2 = gridspec.GridSpec(2, 3, wspace=0.2, hspace=0.05)  ##9 Main Subplots, each divided in 2
    fig2.text(0.5, 0.06, r'Distance [m]', ha='center', fontsize=22)
    fig2.text(0.05, 0.5, r'z [m]', va='center', rotation='vertical', fontsize=22)
    ###Parameters of figure 3 (only npo slid case)
    fig3 = plt.figure(figsize=(16, 8))
    outer3 = gridspec.GridSpec(2, 3, wspace=0.2, hspace=0.05)  ##9 Main Subplots, each divided in 2
    fig3.text(0.5, 0.06, r'Distance [m]', ha='center', fontsize=22)
    fig3.text(0.05, 0.5, r'z [m]', va='center', rotation='vertical', fontsize=22)
    ################################################
    ########## START LOOP ON TRANSECTS##############
    ################################################
    for k, Transect in enumerate(Transect_Names):
        ### Get proper main subplot (one per transect) -> this must be done for fig3 and fig4
        axes1 = np.empty(shape=(1, 3), dtype=object)  ##This will be for subsubplots
        inner1 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer1[k], hspace=0, wspace=0)
        axes2 = np.empty(shape=(1, 3), dtype=object)  ##This will be for subsubplots
        inner2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer2[k], hspace=0, wspace=0)
        axes3 = np.empty(shape=(1, 3), dtype=object)  ##This will be for subsubplots
        inner3 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer3[k], hspace=0, wspace=0)
        ###Get the bed and top DEM of the domain contour
        file_name_bed = './RunTransient_Closure/DEMs/DEM_BedRefined_{}.dat'.format(Transect)
        # Load DEMs of domain
        Col_Names = ['x', 'y']
        DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names,delim_whitespace=True)
        ### Remove first line which contains number of line
        DEM_bed.drop(index=DEM_bed.index[0], axis=0, inplace=True)
        ###Get ybed at x=25 for placing the cross
        df_ybed = DEM_bed.iloc[(DEM_bed['x'] - 25).abs().argsort()[:2]]
        ybedmiddle = Interpolate_on_bed(DEM_bed, 25)

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~~~~~~         GET TUNNEL CONTOUR     ~~~~~~~~###
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###Open the elmer output file
        filename_output_slid = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_Soupirail_{0}_Ovo_Slid_60StepPerY_.dat'.format(Transect)
        filename_output_noslid = 'RunTransient_Closure/RunTransient_{0}_Ovoide/ScalarOutput/TunnelOutput_Soupirail_{0}_Ovo_NoSlid_60StepPerY_.dat'.format(Transect)
        #####Load the Elmer output file
        df0_slid = pd.read_csv(pathroot_mycode.joinpath(filename_output_slid), names=Col_Names_Transient,delim_whitespace=True)
        df0_noslid = pd.read_csv(pathroot_mycode.joinpath(filename_output_noslid), names=Col_Names_Transient,delim_whitespace=True)
        ###################################################
        ########## START LOOP ON OUTPUT TIME ##############
        ###################################################
        for j,(Time,Time_Name) in enumerate(zip(Times,Time_Names)):
            #################################################################################
            ###Get Elmer output corresponding to considered time for the case with sliding###
            df_slid = df0_slid[df0_slid['Tsp']==Time]
            ### get (closed) contour of tunnel for that time
            xt_init_slid = df_slid['X'].values
            yt_init_slid = df_slid['Y'].values
            ###close de contour by re-adding first point
            xt_init_slid = np.append(xt_init_slid, xt_init_slid[0])
            yt_init_slid = np.append(yt_init_slid, yt_init_slid[0])
            ##Sort coordinates of contour clockwise
            xy_slid = np.array([xt_init_slid, yt_init_slid])
            xy_slid = np.transpose(xy_slid)
            df_xy_sorted_slid = pd.DataFrame({'X': xy_slid[:, 0], 'Y': xy_slid[:, 1]})
            ##############################################
            ###Same thing for the case without sliding####
            df_noslid = df0_noslid[df0_noslid['Tsp']==Time]
            ### get (closed) contour of tunnel for that time
            xt_init_noslid = df_noslid['X'].values
            yt_init_noslid = df_noslid['Y'].values
            ###close de contour by re-adding first point
            xt_init_noslid = np.append(xt_init_noslid, xt_init_noslid[0])
            yt_init_noslid = np.append(yt_init_noslid, yt_init_noslid[0])
            ##Sort coordinates of contour clockwise
            xy_noslid = np.array([xt_init_noslid, yt_init_noslid])
            xy_noslid = np.transpose(xy_noslid)
            df_xy_sorted_noslid = pd.DataFrame({'X': xy_noslid[:, 0], 'Y': xy_noslid[:, 1]})
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###~~~         GET PROPER SUB-SUBPLOT AND PLOT          ~~###
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ###Get proper sub-subplots
            ax1 = plt.Subplot(fig1, inner1[j])
            ax2 = plt.Subplot(fig2, inner2[j])
            ax3 = plt.Subplot(fig3, inner3[j])
            ###Name of subsubplots
            for (fig,inner,axes) in zip([fig1,fig2,fig3],[inner1,inner2,inner3],[axes1, axes2, axes3]):
                ###Get proper sub-subplots
                ax = plt.Subplot(fig, inner[j])
                ###Set title of sub-subplot
                ax.set_title('{}'.format(Time_Name), fontsize=16)
                ###Only show axis that are needed
                if k <3:
                    ax.xaxis.set_visible(False) ###xaxis only for bottom row
                if j > 0:
                    ax.yaxis.set_visible(False) ###no yaxis between subsubplots
                ###Set the axe limits and ticks
                xmin = 23.7
                xmax = 26.3
                ###ymin is relative to bed altitude
                df_ybed = DEM_bed.iloc[(DEM_bed['x'] - 25).abs().argsort()[:2]]
                ymin = ybedmiddle - 1.5
                ymax = ybedmiddle + 2.5
                ax.set_xticks([24, 25, 26])
                ax.set_yticks(list(range(2770, 2820, 2)))
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.tick_params(labelsize=16)  # fontsize of the tick labels
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                ax.set_aspect('equal')
                # Add ghost axes and titles on columns
                ax_transect = fig.add_subplot(inner[:])
                ax_transect.axis('off')
                ax_transect.set_title('Transect {}'.format(k+1), fontsize=20, weight='bold')
                ###store current subsubplot
                axes[0, j] = ax
                ###Plot the stuff (both case on some fig)###
                ###Brown color for bed
                ax.fill_between(np.array(DEM_bed['x'].values, dtype=float),np.array(DEM_bed['y'].values, dtype=float), 2750, color='tan')
                ###Icy color for ice
                ax.fill_between(np.array(DEM_bed['x'].values, dtype=float),np.array(DEM_bed['y'].values, dtype=float), 2830, color='lightblue')
                if fig == fig1: ##on fig1 both slid and no slid case
                    ###Plot contour of soupirail and fill in in white for the slid case
                    ax.add_patch(Patch.Polygon(np.transpose(np.array([df_xy_sorted_slid['X'].values,df_xy_sorted_slid['Y'].values])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                    ###print contour in red for no slid cases
                    ax.plot(df_xy_sorted_noslid['X'].values, df_xy_sorted_noslid['Y'].values, color='r', linewidth=2)
                elif fig == fig2: ##on fig2 only slid case
                    ax.add_patch(Patch.Polygon(np.transpose(np.array([df_xy_sorted_slid['X'].values,df_xy_sorted_slid['Y'].values])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                elif fig == fig3: ##on fig3 only no slid case
                    ax.add_patch(Patch.Polygon(np.transpose(np.array([df_xy_sorted_noslid['X'].values,df_xy_sorted_noslid['Y'].values])), closed=True, fill=True, facecolor='w', edgecolor='k'))
                ###print contour of domain
                ax.plot(DEM_bed['x'].values, DEM_bed['y'].values, color='k', linewidth=3)
                ###add current subsubplot to main figure
                fig.add_subplot(ax)
    plt.show()
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###~~~         SAVE THE PLOTS (ONE PER SHAPE/WIDTH)          ~~###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ##Save figure with full transect profile and initial tunnel shape
    name_path_fig = '/home/brondexj/BETTIK/GrandeMotteTignes/RunTransient_Closure/Postprocessing/Figures/.'
    path_output_fig = Path(name_path_fig)
    for (fig, fig_name) in zip([fig1, fig2, fig3], ['Soupirail_SlidAndNoSlid', 'Soupirail_SlidOnly', 'Soupirail_NoSlidOnly']):
        name_output_fig = 'DeformationAt6_9_12months_AllTransects_{0}_Ovoide'.format(fig_name)
        fig.savefig(path_output_fig.joinpath(name_output_fig))






