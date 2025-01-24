"""
@author: jbrondex

Description:
------------
This file aims at evaluating initial cavity volume

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
from scipy.integrate import simpson

from matplotlib import ticker



if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###Open Data corresponding to DEM bed and bottom surf
    Pathroot = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_zbed = 'MNTBedInf2014.dat'
    Filename_zb = 'MNTBedSup2014.dat'
    Col_Names = ['X', 'Y', 'Z']
    ###load files as dataframes
    Df_zbed = pd.read_csv(Pathroot.joinpath(Filename_zbed), names=Col_Names, delim_whitespace=True)
    Df_zb = pd.read_csv(Pathroot.joinpath(Filename_zb), names=Col_Names, delim_whitespace=True)
    ###find min/max x and min/max y
    df_CavityHeight = Df_zb.iloc[:, :2].copy()
    df_CavityHeight['CavityHeight'] = Df_zb['Z']-Df_zbed['Z']
    ###Replace all negative values by zero
    df_CavityHeight.loc[df_CavityHeight['CavityHeight'] < 0, ('CavityHeight')] = 0
    ###Integrate cavity height over surface:
    # Reshape Z into a 2D grid
    x_unique = np.sort(df_CavityHeight['X'].unique())
    y_unique = np.sort(df_CavityHeight['Y'].unique())
    z_grid = df_CavityHeight.pivot(index='Y', columns='X', values='CavityHeight').values
    # Integrate over the surface
    dx = np.diff(x_unique)[0]  # Spacing between X points
    dy = np.diff(y_unique)[0]  # Spacing between Y points
    volume = simpson(simpson(z_grid, x_unique, axis=1), y_unique, axis=0)
    volume_trapz = np.trapz(np.trapz(z_grid, x_unique, axis=1), y_unique, axis=0)
    print(f"Total volume Simpson's rule: {volume}")
    print(f"Total volume trapeze method: {volume_trapz}")
