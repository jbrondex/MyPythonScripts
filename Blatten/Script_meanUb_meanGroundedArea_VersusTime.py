"""
@author: jbrondex

Description:
------------
This file aims at plotting evolution of mean ub and grounded area over time for all Cs

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
# Create legend handles manually, for example:
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path

import numpy as np
import re
import os
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
from matplotlib import ticker
from scipy.interpolate import griddata
from matplotlib.ticker import MultipleLocator

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


###############################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    #### SCRIPT PARAMETERS ###
    ###Parameter of Simu
    Cvalues= [0.3, 0.4, 0.5, 0.6, 0.7]
    # Get 5 shades from the "Blues" colormap (light to dark)
    colors = cmx.get_cmap('viridis_r', len(Cvalues))
    blues = [colors(i) for i in range(len(Cvalues))]
    ###Cases for which we want a figure:
    Cases = ['No Contact', 'Contact']
    LineStyles = ['--','-']
    ###PREPARE THE PLOT
    ###~~~~~~~FIG 1 : mean ub for both with and without contact problem ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax.set_ylabel(r'$\langle u_\mathrm{b} \rangle$ [m/day]', fontsize=20)
    ax.set_xlabel(r'Time [days]', fontsize=20)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ax.set_title('Mean basal velocity [m/day]', fontsize=21, weight='bold')
    ###~~~~~~~FIG 2 : grounded fraction for with contact problem only ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    fig2, ax2 = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    ax2.set_ylabel(r'Grounded Area', fontsize=20)
    ax2.set_xlabel(r'Time [days]', fontsize=20)
    ax2.tick_params(labelsize=16)  # fontsize of the tick labels
    ax2.grid(True)
    ax2.grid(alpha=0.5)
    ax2.set_title('Fraction of grounded area', fontsize=21, weight='bold')

    ### START OF SCRIPT###
    ###Loop on each C considered
    for i, (C, blue_col) in enumerate(zip(Cvalues, blues)):
        print('Dealing with C value :', C)
        str_Cvalue = f"{int(C * 10):02d}"
        ###Loop on cases without and with contact problem
        for j, (Case, LiStyle) in enumerate(zip(Cases, LineStyles)):
            ###Path to Data
            if Case == 'No Contact':
                Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/2_RateAndState_NoDamage_NoContact/GlacierOut/.')
                Filename = 'MyBirch_RS_C{}_WithD_NoCont__bed.dat'.format(str_Cvalue)
                Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta','taub','nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi']
            elif Case == 'Contact':
                Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
                Col_Name = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
                Filename = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue)
                # For following restarts: Regex pattern to extract Restart Number and restart position ($ compulsory to avoid the .dat.names)
                pattern = re.compile(r"MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt(\d+)_tsp(\d+)__bed\.dat$".format(str_Cvalue))
            ###Correct time step size for the simu for which I changes to 0.01 day:
            if 'tsp001d' in Filename:
                Timestep_size = 0.01 / 365  ### In years !!
            else:
                Timestep_size = 0.02 / 365  ### In years !!
            Df_wholebed = pd.DataFrame()
            ###load file as dataframe
            Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
            ###For the case with contact, restarts are required : make the full dataframe from all restarts
            if Case == 'Contact':
                ###number of restart required depends on value of C
                if C == 0.7:
                    Filename_Rstt1 = 'MyBirch_RS_C07_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp2901__bed.dat'
                    Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt1], ignore_index=True)
                elif C == 0.6:
                    ###For this case, there is overlap between end of simu R0 and beginning of simu R1
                    ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 2501)
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2501].reset_index(drop=True)
                    Filename_Rstt1 = 'MyBirch_RS_C06_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp2501__bed.dat'
                    Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt1], ignore_index=True)
                    ###Same thing for restart R2 (tsp 1000 of R1 so 3501 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 3501].reset_index(drop=True)
                    Filename_Rstt2 = 'MyBirch_RS_C06_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp1000__bed.dat'
                    Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt2], ignore_index=True)

                elif C == 0.5:
                    ###For this case, there is overlap between end of simu R0 and beginning of simu R1
                    ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 1601)
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1601].reset_index(drop=True)
                    Filename_Rstt1 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp1601__bed.dat'
                    Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt1], ignore_index=True)
                    ###Same thing for restart R2 (tsp 600 of R1 so 2201 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2201].reset_index(drop=True)
                    Filename_Rstt2 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp600__bed.dat'
                    Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt2], ignore_index=True)
                    ###Same thing for restart R3 (tsp 500 of R2 so 2701 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2701].reset_index(drop=True)
                    Filename_Rstt3 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt3_tsp500__bed.dat'
                    Df_wholebed_Rstt3 = pd.read_csv(Pathroot.joinpath(Filename_Rstt3), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt3 = Df_wholebed_Rstt3[Df_wholebed_Rstt3["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt3["Timestep"] = Df_wholebed_Rstt3["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt3["Count"] = Df_wholebed_Rstt3["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt3], ignore_index=True)
                    ###Same thing for restart R4 (tsp 700 of R3 so 3401 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 3401].reset_index(drop=True)
                    Filename_Rstt4 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt4_tsp700__bed.dat'
                    Df_wholebed_Rstt4 = pd.read_csv(Pathroot.joinpath(Filename_Rstt4), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt4 = Df_wholebed_Rstt4[Df_wholebed_Rstt4["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt4["Timestep"] = Df_wholebed_Rstt4["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt4["Count"] = Df_wholebed_Rstt4["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt4], ignore_index=True)
                elif C == 0.4:
                    ###For this case, there is overlap between end of simu R0 and beginning of simu R1
                    ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 1101)
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1101].reset_index(drop=True)
                    Filename_Rstt1 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp1101__bed.dat'
                    Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt1], ignore_index=True)
                    ###Same thing for restart R2 (tsp 400 of R1 so 1501 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1501].reset_index(drop=True)
                    Filename_Rstt2 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp400__bed.dat'
                    Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt2], ignore_index=True)
                    ###Same thing for restart R3 (tsp 600 of R2 so 2101 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2101].reset_index(drop=True)
                    Filename_Rstt3 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt3_tsp600__bed.dat'
                    Df_wholebed_Rstt3 = pd.read_csv(Pathroot.joinpath(Filename_Rstt3), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt3 = Df_wholebed_Rstt3[Df_wholebed_Rstt3["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt3["Timestep"] = Df_wholebed_Rstt3["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt3["Count"] = Df_wholebed_Rstt3["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt3], ignore_index=True)
                    ###Same thing for restart R4 (tsp 500 of R3 so 2601 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2601].reset_index(drop=True)
                    Filename_Rstt4 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt4_tsp500__bed.dat'
                    Df_wholebed_Rstt4 = pd.read_csv(Pathroot.joinpath(Filename_Rstt4), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt4 = Df_wholebed_Rstt4[Df_wholebed_Rstt4["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt4["Timestep"] = Df_wholebed_Rstt4["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt4["Count"] = Df_wholebed_Rstt4["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt4], ignore_index=True)
                elif C == 0.3:
                    ###For this case, there is overlap between end of simu R0 and beginning of simu R1
                    ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 701)
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 701].reset_index(drop=True)
                    Filename_Rstt1 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp701__bed.dat'
                    Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt1], ignore_index=True)
                    ###Same thing for restart R2 (tsp 300 of R1 so 1001 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1001].reset_index(drop=True)
                    Filename_Rstt2 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp300__bed.dat'
                    Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt2], ignore_index=True)
                    ###Same thing for restart R3 (tsp 500 of R2 so 1501 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1501].reset_index(drop=True)
                    Filename_Rstt3 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt3_tsp500__bed.dat'
                    Df_wholebed_Rstt3 = pd.read_csv(Pathroot.joinpath(Filename_Rstt3), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt3 = Df_wholebed_Rstt3[Df_wholebed_Rstt3["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt3["Timestep"] = Df_wholebed_Rstt3["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt3["Count"] = Df_wholebed_Rstt3["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt3], ignore_index=True)
                    ###Same thing for restart R4 (tsp 400 of R3 so 1901 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 1901].reset_index(drop=True)
                    Filename_Rstt4 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt4_tsp400__bed.dat'
                    Df_wholebed_Rstt4 = pd.read_csv(Pathroot.joinpath(Filename_Rstt4), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt4 = Df_wholebed_Rstt4[Df_wholebed_Rstt4["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt4["Timestep"] = Df_wholebed_Rstt4["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt4["Count"] = Df_wholebed_Rstt4["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt4], ignore_index=True)
                    ###Same thing for restart R5 (tsp 400 of R4 so 2301 from R0) :
                    Df_wholebed = Df_wholebed[Df_wholebed['Timestep'] <= 2301].reset_index(drop=True)
                    Filename_Rstt5 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt5_tsp400__bed.dat'
                    Df_wholebed_Rstt5 = pd.read_csv(Pathroot.joinpath(Filename_Rstt5), names=Col_Name, sep=r'\s+')
                    # 1. Remove rows where Timestep == 1
                    Df_wholebed_Rstt5 = Df_wholebed_Rstt5[Df_wholebed_Rstt5["Timestep"] != 1].copy()
                    # 2. Shift timesteps to continue after the last timestep of Df_wholebed
                    last_timestep = Df_wholebed["Timestep"].max()
                    Df_wholebed_Rstt5["Timestep"] = Df_wholebed_Rstt5["Timestep"] + last_timestep - 1
                    last_count = Df_wholebed["Count"].max()
                    Df_wholebed_Rstt5["Count"] = Df_wholebed_Rstt5["Count"] + last_count - 1
                    # 3. Concatenate back into Df_wholebed
                    Df_wholebed = pd.concat([Df_wholebed, Df_wholebed_Rstt5], ignore_index=True)
            ###Create a column with time (in years !)
            Df_wholebed['Time'] = Timestep_size * Df_wholebed['Timestep']
            ###Calculate ub from velocity and normal vector
            ###Start by calculating un (should be zero but there is some residual)
            Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
            ###Deduce ub from u and un
            Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
            ### Calculate mean ub over whole glacier for each time step
            mean_ub_by_time = Df_wholebed.groupby('Time')['ub'].mean()
            ### Show mean ub in Fig 1
            if LiStyle == '-':
                ax.plot(mean_ub_by_time.index*365, mean_ub_by_time/365, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=2.4)
                ax.legend(fontsize=22)
            else:
                ax.plot(mean_ub_by_time.index * 365, mean_ub_by_time / 365, color=blue_col, linestyle=LiStyle, linewidth=2.4)
            ax.set_ylim(0.0, 40.0)
            ax.set_xlim(0, 24.01)
            # Major ticks at specific values
            ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
            # Minor ticks every 1 unit
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)
            ###Only for the case with problem contact solved, plot evolution of grounded area over time
            if Case == 'Contact':
                ### Calculate total area over the whole glacier for each time step
                Total_Area = Df_wholebed.groupby('Time')['nodearea'].sum()
                ### Grounded area
                Total_GroundedArea = (Df_wholebed[Df_wholebed['GM'] == 1].groupby('Time')['nodearea'].sum())
                ### Fraction of grounded area over total area
                Ratio = Total_GroundedArea/Total_Area
                ### At some precise timestep, ice reground: this seems to be artefact (Free Surface solver missing) so don't plot these points:
                # Copy ratio but mask values where Ratio == 1
                Ratio_masked = Ratio.copy()
                Ratio_masked[Ratio_masked == 1] = np.nan

                ### Plot ratio
                ax2.plot(Total_Area.index*365, Ratio_masked, label= 'C = {}'.format(C), color=blue_col,linestyle=LiStyle, linewidth=2.4)
                ax2.legend(fontsize=22)
                ax2.set_xlim(0, 24.01)
                # Major ticks at specific values
                ax2.set_xticks([0, 4, 8, 12, 16, 20, 24])
                # Minor ticks every 1 unit
                ax2.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax2.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax2.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)




    plt.show()