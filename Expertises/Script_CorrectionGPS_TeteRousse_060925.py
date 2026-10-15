"""
@author: jbrondex

Description:
------------
This file aims at correcting the GPS coordinates measured on Tete Rousse glacier during the visit of 06 September 25.
Error is due to the fact that the wrong plaquette was used for the pivot (plaquette Nord used instead of Pivot Plaquette Sud)

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import matplotlib.dates as mdates


from pathlib import Path

import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
from datetime import datetime, date, timedelta ###To handle dates
import matplotlib.dates as mdates


# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-3, 3))

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
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=20)

###############################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    Pathroot = Path('/home/brondexj/Documents/Expertises/Tete Rousse/Visite06Sept25/.')

    # List of filenames in the desired order
    filename = "TETROUSS06SEPT25_NotCorrected.txt"
    # Define own column names
    column_names = ['Pt Name', 'X', 'Y', 'Z']
    ###Shift of plaquette S (used as ref in the DGPS) relative to plaquette N (where the pivot was placed)
    Xcorr = 0.292
    Ycorr = - 1.196
    Zcorr = - 0.169

    # Read each file while skipping the first row and setting your own column names
    df = pd.read_csv(Pathroot.joinpath(filename), header=None, names=column_names, sep=',')

    df_corrected = df.copy()
    ###Remove the 5 last lines that were not measured on site but already set in the GPS
    df_corrected =df_corrected.iloc[:-6]


    ####Apply the correction
    df_corrected['X'] = df_corrected['X'] - Xcorr
    df_corrected['Y'] = df_corrected['Y'] - Ycorr + 2000000
    df_corrected['Z'] = df_corrected['Z'] - Zcorr


    ###Save the result as a single csv
    Filename_ToSave= 'TETROUSS06SEPT25_Corrected.txt'
    df_corrected.to_csv(Pathroot.joinpath(Filename_ToSave), sep=',', index=False, header=False)



