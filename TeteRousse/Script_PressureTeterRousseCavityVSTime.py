################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path

import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files

from matplotlib import ticker
from scipy import interpolate

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))

###Defined Colors optimized for color-blind people:
RedStar = [162 / 255, 20 / 255, 47 / 255]
RedBackGround = [233 / 255, 201 / 255, 207 / 255]
LightBlueGreyBackGround = [194 / 255, 209 / 255, 219 / 255]
BlueGreyBackGround = [147/ 255, 174 / 255, 191 / 255]
LightBrownBackGround = [223/ 255, 211 / 255, 205 / 255]
VeryLightBrown = [255 / 255, 199 / 255, 127 / 255]
LightBrown = [229 / 255, 143 / 255, 91 / 255]
Brown = [140 / 255, 88 / 255, 56 / 255]
DarkBrown = [64 / 255, 40 / 255, 25 / 255]
VeryDarkBrown = [13 / 255, 8 / 255, 5 / 255]

####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN PART OF THE CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

###Fig 1 is Water level as a funciton of time
fig1 = plt.figure(1, figsize=(40, 20))
plt.ylabel(r'Water Level (m)', fontsize=34)
plt.xlabel(r'Day of Simu', fontsize=34)
plt.xlim([0, 65])
plt.tick_params(labelsize=24)  # fontsize of the tick labels
plt.grid(True)

###############################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    Pathroot = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename = 'Evol_Niveau2010-2013-moyen.dat'
    Col_Names = ['Day', 'WaterLevel']
    Df_WaterLevel = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Names, delim_whitespace=True, skiprows=1)
    ### Renumber days of simu: day 1 of Simu is day -315 of the Water Level File
    Df_WaterLevel['Day']=Df_WaterLevel['Day']+316
    ### Simu start day one and give one output every day:
    DayOfSimu = np.linspace(1,60,60)
    ### Values of Water Level on Simu days are interpolated from water level measured every now and then (Data file)
    WaterLevel_Interp=np.interp(DayOfSimu,Df_WaterLevel['Day'],Df_WaterLevel['WaterLevel']) ##Linear Interpolation as in Elmer
    ###ColorMap for correspondance between Water Pressure and stress
    fig1 = plt.figure(1)
    plt.plot(DayOfSimu,WaterLevel_Interp,color='k',linestyle='-', linewidth=1.5)

    cm = plt.get_cmap('RdBu_r')
    cNorm = colors.Normalize(vmin=np.min(DayOfSimu), vmax=np.max(DayOfSimu))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    plt.scatter(1, WaterLevel_Interp[0], facecolor=scalarMap.to_rgba(1), edgecolors = 'k', s=600, zorder=3)
    for j, day in enumerate(DayOfSimu):
        ### for pumping step profile every 5 days
        if not day % 5 == 0:
            continue
        ###Scatter Plot of Water Level for considered day
        if not day == 60:
            plt.scatter(day, WaterLevel_Interp[j], facecolor=scalarMap.to_rgba(day), edgecolors = 'k', s=600, zorder=3)
        else:
            plt.scatter(day, WaterLevel_Interp[j], facecolor=scalarMap.to_rgba(day), marker='*', edgecolors = 'k', s=1200, zorder=3)
    # plt.axhline(y=SigmathIce, color='k', linestyle=LineStyle, linewidth=2)
    plt.show()