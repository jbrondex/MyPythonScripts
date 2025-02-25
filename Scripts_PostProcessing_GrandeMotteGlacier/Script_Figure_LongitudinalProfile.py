"""
@author: jbrondex

Description:
------------
This scripts is used to produce figure showing the profile of the tunnel along the line
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
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

###Fig 1 is Sigmath as function of RelDens for different parameterization
fig1 = plt.figure(1, figsize=(40, 20))
plt.ylabel(r'z [m]', fontsize=34)
plt.xlabel(r'Distance [m]', fontsize=34)

plt.tick_params(labelsize=24)  # fontsize of the tick labels
plt.grid(True)


if __name__ == "__main__":
#### USER DEFINED PARAMETERS :####
### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/GrandeMotteTignes/.')
    Floor_altitudes = [2803.5, 2803, 2802.5, 2802, 2801.5, 2801.12]
    Transect_xPositions = [47, 97, 147, 197, 252, 291]
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~         CONSTRUCT DOMAIN CONTOUR          ~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###Get the bed and top DEM of the domain contour
    file_name_top = './Maillages/DEMs/DEM_Top_Longitudinal.dat'
    file_name_bed = './Maillages/DEMs/DEM_Bed_Longitudinal.dat'
    # Load DEMs of domain
    Col_Names = ['x', 'y']
    ###DEM top
    DEM_top = pd.read_csv(pathroot_mycode.joinpath(file_name_top), names=Col_Names, delim_whitespace=True)
    Npttop = len(DEM_top.index)
    ###DEM bed
    DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
    ###Missing bed data
    first_row = DEM_bed.iloc[0]
    DEM_bed_missing = pd.concat([first_row.to_frame().T, first_row.to_frame().T], ignore_index=True)
    DEM_bed_missing.iloc[0]['x']=0

    # ###Reverse DEM_bed dataframe so that nodes are consecutive
    # DEM_bed = DEM_bed[::-1]
    # DEM_bed = DEM_bed.reset_index(drop=True)
    # Nptbed = len(DEM_bed.index)


    fig1 = plt.figure(1)
    plt.xlim([DEM_top['x'].min(), DEM_top['x'].max()])
    plt.ylim([2770, 2832])
    ###Brown color for bed
    plt.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2750,color='tan')
    plt.fill_between(np.array(DEM_bed_missing['x'].values, dtype=float), np.array(DEM_bed_missing['y'].values, dtype=float), 2750,color='tan')
    ###Icy color for ice
    plt.fill_between(np.array(DEM_top['x'].values, dtype=float), np.array(DEM_top['y'].values, dtype=float), 2799,color='lightblue')
    plt.fill_between(np.array(DEM_bed['x'].values, dtype=float), np.array(DEM_bed['y'].values, dtype=float), 2801,color='lightblue')
    plt.fill_between(np.array(DEM_bed_missing['x'].values, dtype=float), np.array(DEM_bed_missing['y'].values, dtype=float), 2801,color='lightblue')
    ###White color for tunnel
    plt.fill_between([0, DEM_top['x'].max()], [2804, -0.01*DEM_top['x'].max()+2804 ],  [2804+2.5, -0.01*DEM_top['x'].max()+2804+2.5 ],color='grey')
    # plt.fill_between([0, DEM_top['x'].max()], [2804 + 2.5, -0.01 * DEM_top['x'].max() + 2804 + 2.5], [2786, -0.01 * DEM_top['x'].max() + 2786], color='gainsboro')
    ###Black line for bed
    plt.plot(DEM_bed['x'].values, DEM_bed['y'].values, color='k', linewidth=3)
    plt.plot(DEM_bed_missing['x'].values, DEM_bed_missing['y'].values, color='k', linewidth=3, linestyle='--')
    ###Thin black line for bed
    plt.plot(DEM_top['x'].values, DEM_top['y'].values, color='k', linewidth=1.5)
    ###Thin black line for tunnel
    plt.plot([0, DEM_top['x'].max()], [2804, -0.01*DEM_top['x'].max()+2804 ], color='k', linewidth=1.5)
    plt.plot([0, DEM_top['x'].max()], [2804+2.5, -0.01*DEM_top['x'].max()+2804+2.5 ], color='k', linewidth=1.5)
    ###Thin black line for chenal incision partiel
    plt.plot([0, DEM_top['x'].max()], [2786, -0.01 * DEM_top['x'].max() + 2786], color='k', linewidth=3.2, linestyle=':')
    ###Thick red vertical lines for Transect
    for i in range(len(Transect_xPositions)):
        plt.axvline(x=Transect_xPositions[i], color='red', linestyle='-', linewidth=3)
        plt.text(Transect_xPositions[i]-0.1, 2831.75, 'T{}'.format(i+1), rotation=90, ha='right', va='top', color='red', fontsize=18)
    ##Black vertical line for the transition in between the two bilinear section
    plt.axvline(x=249, color='blue', linestyle=':', linewidth=2)

    plt.show()