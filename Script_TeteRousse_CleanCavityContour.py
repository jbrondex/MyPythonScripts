"""
@author: jbrondex

Description:
------------
This file aims at cleaning cavity contour

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
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###Open Data corresponding to grounded mask
    Pathroot_GM = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Filename_GM = 'GroundedMaskInit.dat'
    Col_Names_GM = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'GM']
    ###load file as dataframe
    Df_GM = pd.read_csv(Pathroot_GM.joinpath(Filename_GM), names=Col_Names_GM, delim_whitespace=True)
    Df_GM.drop_duplicates(inplace=True)
    ###Store the node number of nodes corresponding to GL
    Df_GL=Df_GM[Df_GM['GM']==0]
    ###List of nodes to discard:
    NodesToDiscard=[1506, 202, 978, 1150, 2436, 920, 2274, 2480, 1157, 2022,
                    2957, 2698, 896, 2092, 2234, 196, 2714, 2229, 2421, 908,
                    2830, 2535, 2005, 2188, 927, 2222, 191, 1024, 2254, 1037,
                    1594, 1069, 482, 1991, 138, 2258, 2972, 1232, 2122, 133,
                    2230, 2543, 2141, 2590, 1038, 2471, 2595, 1685, 131, 895,
                    434, 1867, 2202, 2433, 545, 2041, 2252]
    # Create the new DataFrame by filtering out nodes in the NodesToDiscard list
    Df_GL_cleaned = Df_GL[~Df_GL['NodeNumber'].isin(NodesToDiscard)]
    Filename_GL= 'GroundedLine_CLEANED.dat'
    Df_GL_cleaned.to_csv(Pathroot_GM.joinpath(Filename_GL), sep='\t', index=False, header=False)
    # Save the header to a separate file.dat.names
    with open(Pathroot_GM.joinpath('{}.names'.format(Filename_GL)), 'w') as f:
        for i, col in enumerate(Df_GL_cleaned.columns):
            f.write(f"{i+1}: {col}\n")

    ###Sort coordinates of cavity contour clockwise
    xy = np.array([Df_GL_cleaned['X'].values, Df_GL_cleaned['Y'].values])
    xy = np.transpose(xy)
    xy_sorted = sort_clockwise(xy)



    ###Fig 1 is the map of Teterousse, pimp style
    fig1 = plt.figure(1, figsize=(40, 20))
    plt.xlabel(r'X [km]', fontsize=21)
    plt.ylabel(r'Y [km]', fontsize=21)
    plt.tick_params(labelsize=20)  # fontsize of the tick labels
    plt.grid(True)
    plt.grid(alpha=0.5)
    plt.xlim([947950, 948050])
    plt.ylim([2105000,2105120])
    ###Force x and y axis to same scale
    ax = plt.gca()###Plot cavity contour
    ax.set_aspect('equal', adjustable='box')# plt.plot(xy_sorted[:, 0]/1000,xy_sorted[:, 1]/1000,color='dimgrey',linestyle='-',linewidth=1.5)
    plt.scatter(Df_GL['X'].values,Df_GL['Y'].values,color='dimgrey',marker='.',linewidths=0.4)
    # Add labels to each scatter point
    for i, label in enumerate(Df_GL['NodeNumber']):  # Iterate through each point and label    plt.plot(xy_sorted[:, 0],xy_sorted[:, 1],color='k',linestyle='-',linewidth=1)
        plt.text(
            Df_GL['X'].iloc[i],
            Df_GL['Y'].iloc[i],
            label,
            fontsize=9,  # Adjust label font size
            ha='right',  # Horizontal alignment (aligns text to the right of the point)
            va='bottom'
        )
    plt.plot(xy_sorted[:, 0],xy_sorted[:, 1],color='k',linestyle='-',linewidth=1)
    plt.show()