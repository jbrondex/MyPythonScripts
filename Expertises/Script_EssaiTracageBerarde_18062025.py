"""
@author: jbrondex

Description:
------------
This file aims at plotting the measured concentration (mV as calibration was not done) as a function of time
for when we did the rhodamine test during 2025 lake drainage on Bonne Pierre Glacier. Field campaign was on
18/06/2025. They were two fluorimeter.

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
    RefStartDate = datetime(2011, 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    SimuStartDay = -315 ###Number of the day corresponding to the start of the simu pumping 2010
    Pathroot = Path('/home/brondexj/Documents/Expertises/Bérarde/VisiteEssaiTracage_18062025/.')

    # List of filenames in the desired order
    filenames = [
        "BonnePierre_Rhodamine_2025_06_18_1.csv",
        "BonnePierre_Rhodamine_2025_06_18_2.csv",
        "BonnePierre_Rhodamine_2025_06_18_3.csv"
    ]
    # Define own column names
    column_names = ['Day Time', 'Id', 'Value_RG', 'Value_RD', 'Temp_RG', 'Temp_RD', 'None']

    # Define the remarkable time (UTC) for the vertical line and grey shaded
    injection_time = pd.to_datetime('18/06/2025 09:41:00', format='%d/%m/%Y %H:%M:%S')
    BothTabletteRG = pd.to_datetime('18/06/2025 13:28:00', format='%d/%m/%Y %H:%M:%S')
    BothTabletteSamePositionRG = pd.to_datetime('18/06/2025 13:41:00', format='%d/%m/%Y %H:%M:%S')
    BothTabletteRD = pd.to_datetime('18/06/2025 13:57:00', format='%d/%m/%Y %H:%M:%S')
    BothTabletteRD_Julien = pd.to_datetime('18/06/2025 13:58:30', format='%d/%m/%Y %H:%M:%S')
    BothTabletteRD_Gag = pd.to_datetime('18/06/2025 14:01:32', format='%d/%m/%Y %H:%M:%S')
    BothTabletteRD_GagTshirt = pd.to_datetime('18/06/2025 14:02:00', format='%d/%m/%Y %H:%M:%S')

    # Read each file while skipping the first row and setting your own column names
    dfs = [pd.read_csv(Pathroot.joinpath(f), skiprows=1, header=None, names=column_names, sep=';') for f in filenames]
    # Concatenate all dataframes
    df_combined = pd.concat(dfs, ignore_index=True)
    ### Keep only columns that are of interest
    df_measure = df_combined.iloc[:, :4]
    ### Replace 'Non Numérique' by nan
    df_measure.replace('Non Numérique', np.nan, inplace=True)

    # Convert time to datetime.time or datetime.datetime
    df_measure['Day Time'] = pd.to_datetime(df_measure['Day Time'], format='%d/%m/%Y %H:%M:%S')
    # The clock of the tablette was wrong of 3 minutes UTC (tablette was showing 10h02 when it was 9h59 UTC). Correct that :
    df_measure['Day Time'] = df_measure['Day Time'] - pd.Timedelta(minutes=3)
    # Reorder ID column so that it does not start from 0 for each concatenated file
    df_measure['Id'] = range(len(df_measure))

    ###Save the result as a single csv
    Filename_ToSave= 'Tracage_BonnePierre_20250618_Corrected.csv'
    # Define the exact column order and names
    df_to_save = df_measure.rename(columns={
        "Day Time": "Date/Heure",
        "Id": "Point Mesure",
        "Value_RG": "Fluo 1 (valeur D) [mv]",
        "Value_RD": "Fluo 2 (valeur G) [mv]"
    })
    Col_ToSave = ["Date/Heure", "Point Mesure", "Fluo 1 (valeur D) [mv]", ]
    df_to_save.to_csv(Pathroot.joinpath(Filename_ToSave), sep=';', index=False)

    ###ColorMap for correspondance between Water Pressure and stress
    fig1 = plt.figure(1, figsize=(40, 20))
    plt.tick_params(labelsize=18)  # fontsize of the tick labels
    # Define start and end of x-axis manually
    start_time = pd.to_datetime('18/06/2025 09:30:00', format='%d/%m/%Y %H:%M:%S')
    end_time = pd.to_datetime('18/06/2025 14:15:00', format='%d/%m/%Y %H:%M:%S')
    ###Plot the measure RG in red
    plt.plot(df_measure['Day Time'], df_measure['Value_RG'], 'r', linewidth = 1.2, label='Fluo 1')
    plt.plot(df_measure['Day Time'], df_measure['Value_RD'], 'b', linewidth = 1.2, label='Fluo 2')

    # Format x-axis ticks every 15 minutes
    ax = plt.gca()
    ax.set_xlim(start_time, end_time)
    plt.ylim([0, 600])

    # Create ticks every 15 min from 10:00 to last data time
    tick_locs = pd.date_range(start=start_time, end=end_time, freq='15min')
    ax.set_xticks(tick_locs)
    ax.set_xticklabels([t.strftime('%H:%M') for t in tick_locs])

    plt.xticks(rotation=45)
    # fig1.autofmt_xdate()

    ###Shade areas for the different relevant stuff
    # Shade between when fluorimeters are in RG and RD
    plt.axvspan(injection_time, BothTabletteRG, color='turquoise', alpha=0.3, label = 'Fluo 1 RG, Fluo 2 RD')
    # Both Fluo RG but not exact same position
    plt.axvspan(BothTabletteRG, BothTabletteSamePositionRG, color='tomato', alpha=0.3, label = 'Fluo 1 RG, Fluo 2 RG, décalés')
    # Both Fluo RG and exact same position
    plt.axvspan(BothTabletteSamePositionRG, BothTabletteRD, color='r', alpha=0.3, label = 'Fluo 1 RG, Fluo 2 RG, même position')
    # Both Fluo RD and exact same position
    plt.axvspan(BothTabletteRD,end_time, color='grey', alpha=0.3, label = 'Fluo 1 RD, Fluo 2 RD')

    # Plot the vertical line for injection time
    plt.axvline(injection_time, color='k', linestyle='-', linewidth=2)
    # Annotate it with vertical text
    plt.text(injection_time,
             450,  # y_position
             'Injection lac',
             color='k',
             fontsize=16,
             verticalalignment='top',
             horizontalalignment='right',
             rotation=90)

    # Plot a star when Julien washes is hand
    plt.plot(BothTabletteRD_Julien, 580, color='k', marker='*', markersize=18)
    # Annotate it with vertical text
    plt.text(BothTabletteRD_Julien- pd.Timedelta(minutes=2),
             665,  # y_position
             'Mains Julien',
             color='k',
             fontsize=16,
             verticalalignment='top',
             horizontalalignment='right',
             rotation = -45)

    # Plot a star when Gag washes is hand
    plt.plot(BothTabletteRD_Gag, 179, color='k', marker='d', markersize=12)

    # Annotate it with vertical text
    plt.text(BothTabletteRD_GagTshirt + pd.Timedelta(minutes=20),
             259,  # y_position
             'Mains Gag',
             color='k',
             fontsize=16,
             verticalalignment='top',
             horizontalalignment='right',
             rotation = 45)

    # Plot a star when Gag washes is Tshirt
    plt.plot(BothTabletteRD_GagTshirt, 590, color='k', marker='*', markersize=18)
    # Annotate it with vertical text
    plt.text(BothTabletteRD_GagTshirt + pd.Timedelta(minutes=20),
             670,  # y_position
             'Tshirt Gag',
             color='k',
             fontsize=16,
             verticalalignment='top',
             horizontalalignment='right',
             rotation = 45)


    # Final touches
    plt.ylabel(r'Mesure [mV]', fontsize=24)
    plt.xlabel(r'Heure UTC', fontsize=24)
    plt.grid(True)
    plt.title('Lac de Bonne Pierre, Essai de Traçage, 18/06/2025', fontsize=26)
    plt.legend(loc='upper left', bbox_to_anchor=(0.04, 1))
    plt.show()