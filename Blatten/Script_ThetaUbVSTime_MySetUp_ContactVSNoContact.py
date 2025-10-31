"""
@author: jbrondex

Description:
------------
This file aims at plotting several variable (taub, theta, theta-dagger, ...) at a prescribed point of the bed
versus time of simu for my setup

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
    ### In this script I want to compare with and without contact at given C value
    C = 0.3
    str_Cvalue = f"{int(C * 10):02d}"
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 40 ##in m/d !
    ###Path to Data
    Pathroot1 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/2_RateAndState_NoDamage_NoContact/GlacierOut/.')
    Filename1 = 'MyBirch_RS_C{}_WithD_NoCont__bed.dat'.format(str_Cvalue)
    Pathroot2 = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
 ###   Filename2 = 'MyBirch_RS_C{}_Cont_SScalInt10_tsp001d_OutVTU1_DVT_R0_FullPicard__bed.dat'.format(str_Cvalue)
    Filename2 = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue)
    Col_Name1 =  ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi']
    Col_Name2 =  ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub', 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']
    ###Line styles for the various case
    Pathroots= [Pathroot1, Pathroot2]
    Filenames = [Filename1, Filename2]
    Col_Names = [Col_Name1, Col_Name2]
    LiStyles = ['--', '-']
    Labels = ['No Contact', 'Contact']
    ###Parameter of Simu
    Timestep_size = 0.02/365 ### In years !!
    ###Coordinate of bed point at which we want to plot the variables vs time (this is done in the other python script)
    xplot = [2630580, 2630620.9, 2630653.6]
    yplot = [1139265, 1139221.7, 1139188.3]
    colplot = ['#00FF66','#8D5FD3','#FFD42A']

    ###PREPARE THE PLOT
    ###Fig1 is subplot of various variables
    fig1, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    ###Fig 2 same as Fig 1 but organized differently
    fig2, axes2 = plt.subplots(4, 1, figsize=(12, 8), sharex=True)
    ####Figure 3 : tb/N=f(ub/N)
    # Create dedicated figure and axis
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    ###Fig 4 same as Fig 2 but comparing FreeSurface Lateral versus no normal velo for C =0.3 only
    if C ==0.3:
        fig4, axes4 = plt.subplots(4, 1, figsize=(12, 8), sharex=True)

    ### START OF SCRIPT###
    ### To have a legend for the case considered instead of the upper right panel
    handles=[]
    ###Loop on each case (contact vs no contact) considered
    for j,(Pathroot, Filename, Col_Name , Label, LiStyle) in enumerate(zip(Pathroots, Filenames, Col_Names, Labels, LiStyles)):
        ###Correct time step size for the simu for which I changes to 0.01 day:
        if 'tsp001d' in Filename:
            Timestep_size = 0.01/365 ### In years !!
        ###Fill up handles for the common legend regarding C values
        handles.append(Line2D([0], [0], color='black', linestyle=LiStyle, linewidth=3, label=Label))
        Df_wholebed = pd.DataFrame()
        ###load file as dataframe
        Df_wholebed = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Name, sep=r'\s+')
        #####For the with contact case, several restart have been required: build a single dataframe from all restart
        if Pathroot == Pathroot2: ##Only for the with contact case
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
        # We look only the 40 first days : Keep only rows where Time <= 40/365 and reset index
        Df_wholebed = Df_wholebed[Df_wholebed['Time'] <= 40 / 365].reset_index(drop=True)
        # ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes:
        # if Pathroot == Pathroot2: ##Only for the case with contact
        #     # Force taub = 0 where GM < 0
        #     Df_wholebed.loc[Df_wholebed["GM"] < 0, "taub"] = 0
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
        ###Deduce ub from u and un
        Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed['Sigma_nn'] = (Df_wholebed['fx'] * Df_wholebed['nx'] + Df_wholebed['fy'] * Df_wholebed['ny'] + Df_wholebed['fz'] * Df_wholebed['nz']) / Df_wholebed['nodearea']
        # List of variables to interpolate
        # variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'Chi', 'ub', 'Sigma_nn']]
        if Pathroot == Pathroot1: ##For the no contact case, no GM
            variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn']]
        elif Pathroot == Pathroot2: ##For the with contact case, we also interpolate the GM to correct the diagnosed taub
            variable_columns = [col for col in Df_wholebed.columns if col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn', 'GM']]
        ###Loop on each of the point at which we want a curve
        for i,(xstar,ystar,colstar) in enumerate(zip(xplot,yplot,colplot)):
            # Prepare result container
            results = []
            # Loop over each time step
            for time_val, group in Df_wholebed.groupby('Time'):
                # Extract coordinates of nodes and values
                points = group[['X', 'Y']].values
                # Interpolate each variable at the desired point
                interpolated_values = {}
                for var in variable_columns:
                    values = group[var].values
                    interp_val = griddata(points, values, (xstar, ystar), method='linear')
                    interpolated_values[var] = interp_val
                # Store result
                interpolated_values.update({'Time': time_val, 'X': xstar, 'Y': ystar})
                results.append(interpolated_values)
            # Build final DataFrame
            interpolated_df = pd.DataFrame(results)
            # Sort by time
            interpolated_df = interpolated_df.sort_values(by='Time').reset_index(drop=True)
            # ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes for the with contact case:
            if Pathroot == Pathroot2: ##Only for the case with contact
                # Force taub = 0 where GM < 0
                interpolated_df .loc[interpolated_df ["GM"] < 0, "taub"] = 0

            #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            #### FROM THERE I START THE PLOT ####
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 1 ####  : Subplots of theta, taub/N, ub, Sigma_nn
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes[0, 0]
            ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,
                    linewidth=2.2)
            ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta_dagger']), color=colstar,
                    linestyle=LiStyle, linewidth=0.8)
            ax.set_ylabel(r'$\theta$', fontsize=20)
            ax.set_ylim(0.0, 1.0)
            #### DUMMY plot for legend
            if LiStyle == '-' and i == 0:
                ax.plot(np.nan, np.nan, label=r'$\theta$', color='k', linewidth=2.2, linestyle=LiStyle)
                ax.plot(np.nan, np.nan, label=r'$\theta_{dagger}$', color='k', linewidth=0.8, linestyle=LiStyle)
            ax.legend()
            ax.grid(True)

            ###Subplot 2 = nothing for the moment
            axes[0, 1].axis('off')

            ###Subplot 3 = taub/N
            ax = axes[1, 0]
            ###Plot horizontal line for value of C
            if C < 0.6:
                ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
                # Add text label above the line on the left side
                ax.text(
                    35, C,  # x and y coordinates in data units
                    'C = {}'.format(C),
                    ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                    fontsize=10, color='gray'
                )
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['taub'] / interpolated_df['N'], color=colstar,
                    linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.2, 0.6)
            ax.grid(True)

            ###Subplot 4 = nothing for the moment -> use empty subplot to place legend
            axes[1, 1].axis('off')
            # Place legend in the empty subplot
            axes[1, 1].legend(
                handles=handles,
                loc='center',
                fontsize=14,
                frameon=False  # optional: remove box
            )

            ###Subplot 5 : ub
            ax = axes[2, 0]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_ylim(0.0, velocitymax)
            ax.grid(True)
            ###Because axes[2,1] is off, put xlabel and force ticks here
            ax.tick_params(labelbottom=True)
            ax.set_xlabel(r'Time [days]', fontsize=20)

            ###Subplot 6 : Sigmann
            ax = axes[2, 1]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.set_ylim(-5.0, 5.0)
            ax.grid(True)

            ####Spare Subplot : Chi
            # ax.plot(interpolated_df['Time']*365, interpolated_df['Chi'], color=colstar, linestyle=LiStyle, linewidth=2)
            # ax.axhspan(0, 1.8, color='grey', alpha=0.05)
            # ax.set_ylabel(r'$\chi$', fontsize=20)
            # ax.set_xlabel(r'Time [days]', fontsize=20)
            # ax.set_ylim(-3.5,1.8)
            # ax.grid(True)

            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 2 ####  : Same as Fig 1 but one single column
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes2[0]
            ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,
                    linewidth=2.2)
            ax.set_ylabel(r'$\theta$', fontsize=20)
            ax.set_ylim(0.0, 1.0)
            ###xlim for all plots : only first 20 days
            ax.set_xlim(0, 24.01)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ###Subplot 2 = taub/N
            ax = axes2[1]
            ###Plot horizontal line for value of C
            if C < 0.6:
                ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
                # Add text label above the line on the left side
                ax.text(
                    16, C,  # x and y coordinates in data units
                    'C = {}'.format(C),
                    ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                    fontsize=10, color='gray'
                )
            ####Everywhere taub is zero I replace by NaN to have a gap in the solid line (makes visualisation easier)
            interpolated_df_taubNaN = interpolated_df.copy()
            interpolated_df_taubNaN.loc[interpolated_df_taubNaN["taub"] == 0, "taub"] = np.nan
            ax.plot(interpolated_df_taubNaN['Time'] * 365, interpolated_df_taubNaN['taub'] / interpolated_df_taubNaN['N'], color=colstar,
                    linestyle=LiStyle, linewidth=2)
            ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
            ax.set_ylim(0.2, 0.7)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)


            ###Subplot 3 : ub
            ax = axes2[2]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
            ax.set_ylim(0.0, velocitymax)
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ###Subplot 4 : Sigmann
            ax = axes2[3]
            ax.plot(interpolated_df['Time'] * 365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,
                    linewidth=2)
            ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
            ax.set_xlabel(r'Time [days]', fontsize=20)
            ax.set_ylim(-5.0, 5.0)
            # Major ticks at specific values
            ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
            # Minor ticks every 1 unit
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            # Add grid for both major and minor ticks
            ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
            ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 3 ####
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            # Fill up figure 3 : Plot the ratio ub/N vs taub/N
            ax3.plot(interpolated_df['ub'] / (365 * interpolated_df['N']),
                     interpolated_df['taub'] / interpolated_df['N'], color=colstar, linestyle=LiStyle, linewidth=2)
            # Axis labels with LaTeX-style formatting
            ax3.set_xlabel(r'$\mathrm{u}_\mathrm{b} / N$ [m d$^{-1}$ MPa$^{-1}$]', fontsize=20)
            ax3.set_ylabel(r'$\tau_b / N$', fontsize=20)
            ax3.grid(True)

            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            #### FIGURE 4 ####  : Same as Fig 2 but for the case of free surf lat (C=0.3 only)
            ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            if C == 0.3:
                ###Subplot 1 : Theta
                ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
                ax = axes4[0]
                ###xlim for all plots : only first 20 days
                ax.set_xlim(0, 24.01)
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ax.plot(interpolated_df['Time'] * 365, (1 - interpolated_df['theta']), color=colstar, linestyle=LiStyle,
                        linewidth=2.2)
                ax.set_ylabel(r'$\theta$', fontsize=20)
                ax.set_ylim(0.0, 1.0)



                ###Subplot 2 = taub/N
                ax = axes4[1]
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ###Plot horizontal line for value of C
                if C < 0.6:
                    ax.axhline(y=C, color='gray', linestyle='--', linewidth=1.5)
                    # Add text label above the line on the left side
                    ax.text(
                        21, C,  # x and y coordinates in data units
                        'C = {}'.format(C),
                        ha='left', va='bottom',  # align text to left and bottom  # place text in axes (0–1) coordinates
                        fontsize=14, color='gray'
                    )
                ####Everywhere taub is zero I replace by NaN to have a gap in the solid line (makes visualisation easier)
                interpolated_df_taubNaN = interpolated_df.copy()
                interpolated_df_taubNaN.loc[interpolated_df_taubNaN["taub"] == 0, "taub"] = np.nan
                ax.plot(interpolated_df_taubNaN['Time'] * 365, interpolated_df_taubNaN['taub'] / interpolated_df_taubNaN['N'], color=colstar,
                        linestyle=LiStyle, linewidth=2)
                ax.set_ylabel(r'$\tau_b$ / N', fontsize=20)
                ax.set_ylim(0.2, 1.0)


                ###Subplot 3 : ub
                ax = axes4[2]
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ax.plot(interpolated_df['Time'] * 365, interpolated_df['ub'] / 365, color=colstar, linestyle=LiStyle,
                        linewidth=2)
                ax.set_ylabel(r'$\mathrm{u}_\mathrm{b}$ [m/d]', fontsize=20)
                ax.set_ylim(0.0, velocitymax+10)


                ###Subplot 4 : Sigmann
                ax = axes4[3]
                # Major ticks at specific values
                ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
                # Minor ticks every 1 unit
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                # Add grid for both major and minor ticks
                ax.grid(which="major", color="black", linestyle="-", linewidth=0.6, alpha=0.6)
                ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.4, alpha=0.6)

                ax.plot(interpolated_df['Time'] * 365, interpolated_df['Sigma_nn'], color=colstar, linestyle=LiStyle,
                        linewidth=2)
                ax.set_ylabel(r'$\sigma_\mathrm{nn}$ [MPa]', fontsize=20)
                ax.set_xlabel(r'Time [days]', fontsize=20)
                ax.set_ylim(-7.0, 2.2)


    #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    #### HERE I TAKE CARE OF THE CASE FREE SURF LATERAL FOR C =0.3 and WITH CONTACT PROBLEM ONLY
    if C == 0.3:
        Filename3 = 'MyBirch_RS_C03_Cont_SScalInt10_tsp001d_OutVTU1_DVT_FullPicard_LatFS__bed.dat'
        Timestep_size = 0.01 / 365  ### In years !!
        ###Fill up handles for the common legend regarding C values
        Df_wholebed3 = pd.DataFrame()
        ###load file as dataframe
        Df_wholebed3 = pd.read_csv(Pathroot2.joinpath(Filename3), names=Col_Name2, sep=r'\s+')
        ###Create a column with time (in years !)
        Df_wholebed3['Time'] = Timestep_size * Df_wholebed3['Timestep']
        ###Calculate ub from velocity and normal vector
        ###Start by calculating un (should be zero but there is some residual)
        Df_wholebed3['un'] = Df_wholebed3['u'] * Df_wholebed3['nx'] + Df_wholebed3['v'] * Df_wholebed3['ny'] + \
                             Df_wholebed3['w'] * Df_wholebed3['nz']
        ###Deduce ub from u and un
        Df_wholebed3['ub'] = np.sqrt((Df_wholebed3['u'] - Df_wholebed3['un'] * Df_wholebed3['nx']) ** 2 + (
                Df_wholebed3['v'] - Df_wholebed3['un'] * Df_wholebed3['ny']) ** 2 + (
                                             Df_wholebed3['w'] - Df_wholebed3['un'] * Df_wholebed3['nz']) ** 2)
        ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
        Df_wholebed3['Sigma_nn'] = (Df_wholebed3['fx'] * Df_wholebed3['nx'] + Df_wholebed3['fy'] * Df_wholebed3[
            'ny'] + Df_wholebed3['fz'] * Df_wholebed3['nz']) / Df_wholebed3['nodearea']
        # List of variables to interpolate
        variable_columns = [col for col in Df_wholebed3.columns if
                            col in ['H', 'N', 'theta_dagger', 'theta', 'taub', 'fz', 'ub', 'Sigma_nn', 'GM']]
        ###Loop on each of the point at which we want a curve
        for i, (xstar, ystar, colstar) in enumerate(zip(xplot, yplot, colplot)):
            # Prepare result container
            results3 = []
            # Loop over each time step
            for time_val, group in Df_wholebed3.groupby('Time'):
                # Extract coordinates of nodes and values
                points = group[['X', 'Y']].values
                # Interpolate each variable at the desired point
                interpolated_values3 = {}
                for var in variable_columns:
                    values = group[var].values
                    interp_val = griddata(points, values, (xstar, ystar), method='linear')
                    interpolated_values3[var] = interp_val
                # Store result
                interpolated_values3.update({'Time': time_val, 'X': xstar, 'Y': ystar})
                results3.append(interpolated_values3)
            # Build final DataFrame
            interpolated_df3 = pd.DataFrame(results3)
            # Sort by time
            interpolated_df3 = interpolated_df3.sort_values(by='Time').reset_index(drop=True)
            # Force taub = 0 where GM < 0
            interpolated_df3 .loc[interpolated_df3 ["GM"] < 0, "taub"] = 0
            ### FILL-UP FIGURE 4 WITH THIS LATERAL FREE SURFACE CASE
            ###Subplot 1 : Theta
            ### BE CAREFUL, WE CHANGE CONVENTION COMPARED TO ELMER : theta = 0 for no cavitation, theta = 1 for full cavitation
            ax = axes4[0]
            ax.plot(interpolated_df3['Time'] * 365, (1 - interpolated_df3['theta']), color=colstar, linestyle=':',linewidth=2.2)
            ###Subplot 2 = taub/N
            ax = axes4[1]
            ####Everywhere taub is zero I replace by NaN to have a gap in the solid line (makes visualisation easier)
            interpolated_df3_taubNaN = interpolated_df3.copy()
            interpolated_df3_taubNaN.loc[interpolated_df3_taubNaN["taub"] == 0, "taub"] = np.nan
            ax.plot(interpolated_df3_taubNaN['Time'] * 365, interpolated_df3_taubNaN['taub'] / interpolated_df3_taubNaN['N'], color=colstar,linestyle=':', linewidth=2)
            ###Subplot 3 : ub
            ax = axes4[2]
            ax.plot(interpolated_df3['Time'] * 365, interpolated_df3['ub'] / 365, color=colstar, linestyle=':',linewidth=2)
            ###Subplot 4 : Sigmann
            ax = axes4[3]
            ax.plot(interpolated_df3['Time'] * 365, interpolated_df3['Sigma_nn'], color=colstar, linestyle=':',linewidth=2)

    ###Now add legend regarding C on empty upper right subplot
    axes[1, 1].legend(handles=handles,loc='center',fontsize=22,frameon=False)

    plt.show()