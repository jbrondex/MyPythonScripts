"""
@author: jbrondex

Description:
------------
This file aims at plotting the evolution of the mean vertical velocity above the cavity between first pumping of 2010 to 2014.

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files

from matplotlib import ticker

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
############################################################
# Define function to check wether node is above cavity #####
############################################################
def IsAboveCavity(x,y):
    xc = 948003
    Rx = 22
    yc = 2105064
    Ry = 54
    metric = ((xc - x)/Rx)**2 + ((yc - y)/Ry)**2
    out = metric <= 1 ###LOGICAL: True if above cavity, False otherwise
    return out

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

fig, axes = plt.subplots(1, 3, figsize=(60, 20), sharey= True)
# pimp style
##First Row
ax = axes[0]
ax.set_ylabel(r'Mean Vertical Velocity (mm/d)', fontsize=34)
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)


ax = axes[1]
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
# ###Plot vertical lines corresponding to min and max rel density
# ax.axvline(x=0.38, color= 'k', linestyle = ':', linewidth=3)
# ax.axvline(x=1, color= 'k', linestyle = ':', linewidth=3)

ax = axes[2]
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
# ###Plot vertical lines corresponding to min and max rel density
# ax.axvline(x=0.38, color= 'k', linestyle = ':', linewidth=3)
# ax.axvline(x=1, color= 'k', linestyle = ':', linewidth=3)


################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    #########################################
    #### FIRST OUTPUT OF THE SIMULATIONS:####
    #########################################
    ### General parameter of simu
    RefStartDate = date(int(2011), 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    StartDayNumber_Pumping2010 = "-315"
    StartYear_Pumping2010 = "2010"
    StartDate_Pumping2010 = RefStartDate + timedelta(days=int(StartDayNumber_Pumping2010))

    ###Load file containing all tested parameter sets for damage model
    Pathroot_ParamSet = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/.')
    Filename_ParamSet = 'Dam_Sigmath_B_Lambdah.IN'
    Col_Names_ParamSet = ['Sigmath', 'B', 'Lambdah']
    Data_ParamSet = pd.read_csv(Pathroot_ParamSet.joinpath(Filename_ParamSet), names=Col_Names_ParamSet, delim_whitespace=True)
    ###Add colums to Data_ParamSet corresponding to parameters converted in file names (convert to string and remove dot)
    Data_ParamSet['Sigmath_name']=Data_ParamSet['Sigmath'].astype(str).replace('\.', '',regex=True)
    Data_ParamSet['B_name']=Data_ParamSet['B'].astype(str).replace('\.', '',regex=True)
    Data_ParamSet['Lambdah_name']=Data_ParamSet['Lambdah'].astype(str).replace('\.', '',regex=True)
    ###List of tested lambdah values
    Lambdah_list = [0.0, 0.5, 1.0]
    ###Corresponding linestyle
    LineStyle_List = ['-', '--', ':']
    ###Create a color map to cover all combinations of (Sigmath, B) for a given value of lamddah
    cm = plt.get_cmap('copper')
    cNorm = colors.Normalize(vmin=0, vmax=len(Data_ParamSet[(Data_ParamSet['Lambdah']==0.0)]))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    Color_List = []
    for i in range(len(Data_ParamSet[(Data_ParamSet['Lambdah']==Lambdah_list[0])])):
        Color_List.append(scalarMap.to_rgba(i))
    ###Alternative choice for colors: a discrete colormap that is colorblind-friendly
    ###In this order: blue/orange/green/pink/brown/purple/grey/red/yellow
    ColorBlindFriendly_List = ['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
    ### Load files containing surface output of the simulations:
    Pathroot_SimuOutput = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    ### Step in the full process from 2010 to 2014
    Step_List = ['Pump2010', 'Refill20102011']
    StepTsp_List =[1, 5] ##Timestep size (in days) corresponding to simulation step (Step_List)
    StepOut_List =[5, 30] ##Output interval (in days) corresponding to simulation step (Step_List)
    ### Where pressure applies: cavity only, restricted to a conduit, everywhere above cold/temperate transition
    Cases = ['PCavityOnly', 'PRestricted', 'PNotRestric']
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'SigmaI', 'Damage', 'Chi']
    Col_Names_NoD = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'SigmaI']  ##For the no damage (ref) case
    ###START LOOP over each pressure scenario (one subplot par scenario)
    for j, case in enumerate(Cases):
        ### Get the corresponding subplot
        ax = axes[j]
        ax.set_title('{}'.format(case), fontsize=34, weight='bold')
        ###For each pressure scenario get the no damage (ref) case
        ###We create a single dataframe for all steps
        Data_Simu_NoD = pd.DataFrame() ##Initialize an empty dataframe
        for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
            filename = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_dm315todm255_Out{}d_.dat'.format(case, Step, str(StepTsp), str(StepOut))
            print('Opening file:', filename)
            Data_Simu_NoD_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names_NoD, delim_whitespace=True)
            ###Drop duplicate lines (bug of SaveLine solver)
            Data_Simu_NoD_tmp.drop_duplicates(inplace=True)
            ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
            if i == 0:
                Data_Simu_NoD_tmp['DayOfSimu'] = Data_Simu_NoD_tmp['DayOfSimu']*StepTsp
            else:
                Data_Simu_NoD_tmp['DayOfSimu'] = np.max(Data_Simu_NoD['DayOfSimu']) + Data_Simu_NoD_tmp['DayOfSimu']*StepTsp
            data = [Data_Simu_NoD,Data_Simu_NoD_tmp]
            Data_Simu_NoD=pd.concat(data, ignore_index=True)
        ###Add a column of boolean to dataset depending on wether node is above cavity or not
        Data_Simu_NoD['IsAboveCavity']=IsAboveCavity(Data_Simu_NoD['X'],Data_Simu_NoD['Y'])
        ###Store Node above cavity in a separated dataframe
        Data_Simu_NoD_AboveCavity=Data_Simu_NoD[Data_Simu_NoD['IsAboveCavity']]
        ###For each timestep calculate mean, min, max of vertical velocity above cavity
        Date = []
        MeanW_NoD = []
        MinW_NoD = []
        MaxW_NoD = []
        ### We want to product a list corresponding to all simulation days at which we have an output
        SimuDays = Data_Simu_NoD['DayOfSimu']
        SimuDays.drop_duplicates(inplace=True)
        for k in range(len(SimuDays)):
            day =  int(SimuDays.iloc[k])
            ###Convert timestep in terms of date
            Date.append(StartDate_Pumping2010 + timedelta(days=day-1))
            ###Calculate mean/min/max in mm/days
            MeanW_NoD.append(np.mean(Data_Simu_NoD_AboveCavity[Data_Simu_NoD_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
            MinW_NoD.append(np.min(Data_Simu_NoD_AboveCavity[Data_Simu_NoD_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
            MaxW_NoD.append(np.max(Data_Simu_NoD_AboveCavity[Data_Simu_NoD_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
        ###Plot MeanW for the case with no damage (ref case) on corresponding subplot
        ax.plot_date(Date, MeanW_NoD, color= 'k', linestyle = '-', linewidth=2, marker='None', xdate=True)
        ax.xaxis.set_major_locator(MonthLocator(interval=2))
        ax.xaxis.set_minor_locator(DayLocator(interval=14))
        ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
        fig.autofmt_xdate()

        ###From there do the loop on each damage parameter set
        for l in range(len(Data_ParamSet)):
            print('Considering Damage parameter set n°',l+1,'/',len(Data_ParamSet))
            ###It turns out that results with lambdah=0.5 are very close to the ones with lambdah=1.0 so skip the latter
            if Data_ParamSet['Lambdah'][l] == Lambdah_list[2]:
                continue
            ###We create a single dataframe for all steps
            Data_Simu_D = pd.DataFrame()  ##Initialize an empty dataframe
            for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
                filename = 'SurfaceOutput_Prono_Dam_B{}_Sigth{}_Lambdah{}_{}_{}_tsp{}d_dm315todm255_Out{}d_.dat'.format(Data_ParamSet['B_name'][l], Data_ParamSet['Sigmath_name'][l], Data_ParamSet['Lambdah_name'][l], case, Step, str(StepTsp), str(StepOut))
                print('Opening file:', filename)
                Data_Simu_D_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names, delim_whitespace=True)
                ###Drop duplicate lines (bug of SaveLine solver
                Data_Simu_D_tmp.drop_duplicates(inplace=True)
                ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
                if i == 0:
                    Data_Simu_D_tmp['DayOfSimu'] = Data_Simu_D_tmp['DayOfSimu'] * StepTsp
                else:
                    Data_Simu_D_tmp['DayOfSimu'] = np.max(Data_Simu_D['DayOfSimu']) + Data_Simu_D_tmp['DayOfSimu'] * StepTsp
                data_1 = [Data_Simu_D, Data_Simu_D_tmp]
                Data_Simu_D = pd.concat(data_1, ignore_index=True)
            ###Add a column of boolean to dataset depending on wether node is above cavity or not
            Data_Simu_D['IsAboveCavity']=IsAboveCavity(Data_Simu_D['X'],Data_Simu_D['Y'])
            ###Store Node above cavity in a separated dataframe
            Data_Simu_D_AboveCavity=Data_Simu_D[Data_Simu_D['IsAboveCavity']]
            ###For each timestep calculate mean, min, max of vertical velocity above cavity
            Date = []
            MeanW_D = []
            MinW_D = []
            MaxW_D = []
            ### We want to product a list corresponding to all simulation days at which we have an output
            SimuDays = Data_Simu_D['DayOfSimu']
            SimuDays.drop_duplicates(inplace=True)
            for k in range(len(SimuDays)):
                day = int(SimuDays.iloc[k])
                ###Convert timestep in terms of date
                Date.append(StartDate_Pumping2010 + timedelta(days=day-1))
                ###Calculate mean/min/max in mm/days
                MeanW_D.append(np.mean(Data_Simu_D_AboveCavity[Data_Simu_D_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
                MinW_D.append(np.min(Data_Simu_D_AboveCavity[Data_Simu_D_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
                MaxW_D.append(np.max(Data_Simu_D_AboveCavity[Data_Simu_D_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
            ###Plot MeanW for the case with no damage (ref case) on corresponding subplot
            if Data_ParamSet['Lambdah'][l] == Lambdah_list[0]:
                LineStyle = LineStyle_List[0]
                Color = ColorBlindFriendly_List[l]
            elif Data_ParamSet['Lambdah'][l] == Lambdah_list[1]:
                LineStyle = LineStyle_List[1]
                Color = ColorBlindFriendly_List[l - 8]
            elif Data_ParamSet['Lambdah'][l] == Lambdah_list[2]:
                LineStyle = LineStyle_List[2]
                Color = ColorBlindFriendly_List[l - 16]
            else:
                print('ERROR: Lambdah_list does not correspond to damage parameter sets input file')
            ax.plot_date(Date, MeanW_D, color= Color, linestyle = LineStyle, linewidth=1.5, marker='None', xdate=True)





    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'RelDensProfile_Forage{}_{}_DInitHL_'.format(ForageNumber,  Case)
    # name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_FullSimu_TestNbIntTsp/'
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()