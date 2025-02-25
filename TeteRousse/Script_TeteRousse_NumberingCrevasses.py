"""
@author: jbrondex

Description:
------------
This file aims at numbering crevasses

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from scipy.interpolate import griddata
from matplotlib.ticker import FormatStrFormatter
from matplotlib.dates import date2num
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files

from matplotlib import ticker

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
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)


if __name__ == "__main__":
    ################################################################
    ####  OPEN DATA CORRESPONDING TO CREVASSES AND GROUNDEDMASK ####
    ################################################################
    Pathroot_Obs = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_Crevasses = 'crevasses_xyz.dat'
    ###load file as dataframe
    Df_Crevasses = pd.read_csv(Pathroot_Obs.joinpath(Filename_Crevasses), names=['X', 'Y', 'Z'], delim_whitespace=True)
    ###Also load glacier contour
    Filename_GlacierContour = 'Contour_TeteRousse2012.dat'
    Df_GlacierContour = pd.read_csv(Pathroot_Obs.joinpath(Filename_GlacierContour), names=['X', 'Y'], delim_whitespace=True)
    xc = Df_GlacierContour['X'].values
    xc=np.append(xc,xc[0]) ##step required to close the contour
    yc = Df_GlacierContour['Y'].values
    yc = np.append(yc, yc[0])  ##step required to close the contour
    ###Open Data corresponding to grounded mask
    Pathroot_GM = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Filename_GM = 'GroundedLine_CLEANED.dat'
    Col_Names_GM = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'GM']
    ###load file as dataframe
    Df_GL = pd.read_csv(Pathroot_GM.joinpath(Filename_GM), names=Col_Names_GM, delim_whitespace=True)
    ###Sort coordinates of cavity contour clockwise
    xy_cavity = np.array([Df_GL['X'].values, Df_GL['Y'].values])
    xy_cavity = np.transpose(xy_cavity)
    xy_cavity_sorted = sort_clockwise(xy_cavity)
    ##close contour
    cavity_contour = np.vstack([xy_cavity_sorted, xy_cavity_sorted[0]])

    #######################################################################################
    #### PROCESSING TO NUMBER THE CREVASSES BASE ON INDEX OF GPS POINTS IN Df_Crevasses ###
    #######################################################################################
    ### Initialise new columns to NaN
    Df_Crevasses['Crevasse Number'] = 0
    Df_Crevasses['IsCircular'] = False
    ###List indexes of Df_Crevasses on a crevasse per crevasse basis
    Df_Crevasses.loc[np.arange(415,423), 'Crevasse Number'] = 1
    Df_Crevasses.loc[np.arange(415,423), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(399,415), 'Crevasse Number'] = 2
    Df_Crevasses.loc[np.arange(399,415), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(359,399), 'Crevasse Number'] = 3
    Df_Crevasses.loc[np.arange(359,399), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(344,359), 'Crevasse Number'] = 4
    Df_Crevasses.loc[np.arange(344,359), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(318,327), 'Crevasse Number'] = 5
    Df_Crevasses.loc[np.arange(318,327), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(327,333), 'Crevasse Number'] = 6
    Df_Crevasses.loc[np.arange(327,333), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(333,344), 'Crevasse Number'] = 7
    Df_Crevasses.loc[np.arange(333,344), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(295,318), 'Crevasse Number'] = 8
    Df_Crevasses.loc[np.arange(295,318), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(277,295), 'Crevasse Number'] = 9
    Df_Crevasses.loc[np.arange(277,295), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(78,88), 'Crevasse Number'] = 10
    Df_Crevasses.loc[np.arange(78,88), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(52,78), 'Crevasse Number'] = 11
    Df_Crevasses.loc[np.arange(52,78), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(23,52), 'Crevasse Number'] = 12
    Df_Crevasses.loc[np.arange(23,52), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(88,97), 'Crevasse Number'] = 13
    Df_Crevasses.loc[np.arange(88,97), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(0,23), 'Crevasse Number'] = 14
    Df_Crevasses.loc[np.arange(0,23), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(97,128), 'Crevasse Number'] = 15
    Df_Crevasses.loc[np.arange(97,128), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(128,152), 'Crevasse Number'] = 16
    Df_Crevasses.loc[np.arange(128,152), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(193,207), 'Crevasse Number'] = 17
    Df_Crevasses.loc[np.arange(193,207), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(183,193), 'Crevasse Number'] = 18
    Df_Crevasses.loc[np.arange(183,193), 'IsCircular'] = True
    Df_Crevasses.loc[np.arange(207,224), 'Crevasse Number'] = 19
    Df_Crevasses.loc[np.arange(224,236), 'Crevasse Number'] = 20
    Df_Crevasses.loc[np.arange(236,246), 'Crevasse Number'] = 21
    Df_Crevasses.loc[np.arange(246,256), 'Crevasse Number'] = 22
    Df_Crevasses.loc[np.arange(256,261), 'Crevasse Number'] = 23
    Df_Crevasses.loc[np.arange(261,277), 'Crevasse Number'] = 24
    ###Remove all lines containing 0 as crevasse number (means that the gps point was not associated to any crevasse)
    Df_Crevasses = Df_Crevasses.loc[Df_Crevasses['Crevasse Number'] != 0]

    ####Save this dataframe for using it elsewhere
    Filename_CrevassesNumbered= 'crevasses_xyz_numbered.dat'
    Df_Crevasses.to_csv(Pathroot_Obs.joinpath(Filename_CrevassesNumbered), sep='\t', index=False, header=False)
    # Save the header to a separate file.dat.names
    with open(Pathroot_Obs.joinpath('{}.names'.format(Filename_CrevassesNumbered)), 'w') as f:
        for i, col in enumerate(Df_Crevasses.columns):
            f.write(f"{i+1}: {col}\n")


    #############################
    #### DO THE PLOT TO CHECK ###
    #############################
    ###Fig 1 is the map of Teterousse, pimp style
    fig1 = plt.figure(1, figsize=(40, 20))
    plt.xlabel(r'X [km]', fontsize=21)
    plt.ylabel(r'Y [km]', fontsize=21)
    plt.tick_params(labelsize=20)  # fontsize of the tick labels
    plt.grid(True)
    plt.grid(alpha=0.5)
    # plt.xlim([947950, 948050])
    # plt.ylim([2105000,2105120])
    ###Force x and y axis to same scale
    ax = plt.gca()
    plt.plot(xc, yc, color='k', linestyle='-', linewidth=1)
    plt.plot(cavity_contour[:, 0], cavity_contour[:, 1], color='k', linestyle='-', linewidth=1)
    ax.set_aspect('equal', adjustable='box')# plt.plot(xy_sorted[:, 0]/1000,xy_sorted[:, 1]/1000,color='dimgrey',linestyle='-',linewidth=1.5)
    plt.scatter(Df_Crevasses['X'].values,Df_Crevasses['Y'].values,color='g',marker='+',linewidths=0.4)
    # Add labels to each scatter point
    for i, label in enumerate(Df_Crevasses.index):  # Iterate through each point and label    plt.plot(xy_sorted[:, 0],xy_sorted[:, 1],color='k',linestyle='-',linewidth=1)
        plt.text(
            Df_Crevasses['X'].iloc[i],
            Df_Crevasses['Y'].iloc[i],
            label,
            fontsize=9,  # Adjust label font size
            ha='right',  # Horizontal alignment (aligns text to the right of the point)
            va='bottom'
        )
    ###plot crevasses as continuous line
    for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(),Df_Crevasses['Crevasse Number'].max()+1):
        ###get proper points
        Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number']==crev_num]
        if Df_plot['IsCircular'].all():##different colors for circular crevasses and other crevasses
            col = 'b'
        else:
            col = 'r'
        plt.plot(Df_plot['X'].values, Df_plot['Y'].values, color=col, linestyle='-', linewidth=1.5)
        # # Add labels corresponding to crevasse number
        # label = int(Df_plot['Crevasse Number'].iloc[0])
        # plt.text(Df_plot['X'].iloc[0],Df_plot['Y'].iloc[0],label,fontsize=9,ha='right', va='bottom')

    plt.show()