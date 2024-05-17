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
from matplotlib.dates import date2num
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

fig, axes = plt.subplots(2, 2, figsize=(60, 20), sharex= True , sharey= True)
# pimp style
##First Row
ax = axes[0,0]
ax.set_ylabel(r'Mean $V_z$ (mm/d)', fontsize=25)
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
ax.xaxis.set_major_locator(MonthLocator(interval=3))
ax.xaxis.set_minor_locator(MonthLocator(interval=1))
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
fig.autofmt_xdate()
###plot one horizontal line for Vz=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)
###Shade periods corresponding to pumping (26/08/2010 to 15/10/10; 28/09/11 to14/10/11; 23/09/12 to 09/10/12)
ax.axvspan(date(int(2010), 8, 26) ,date(int(2010), 10, 15) , alpha=0.3, color='grey')
ax.axvspan(date(int(2011), 9, 28) ,date(int(2011), 10, 14) , alpha=0.3, color='grey')
ax.axvspan(date(int(2012), 9, 23) ,date(int(2012), 10, 9) , alpha=0.3, color='grey')

ax = axes[0,1]
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
ax.xaxis.set_major_locator(MonthLocator(interval=3))
ax.xaxis.set_minor_locator(MonthLocator(interval=1))
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
fig.autofmt_xdate()
###plot one horizontal line for Vz=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)
###Shade periods corresponding to pumping (26/08/2010 to 15/10/10; 28/09/11 to14/10/11; 23/09/12 to 09/10/12)
ax.axvspan(date(int(2010), 8, 26) ,date(int(2010), 10, 15) , alpha=0.3, color='grey')
ax.axvspan(date(int(2011), 9, 28) ,date(int(2011), 10, 14) , alpha=0.3, color='grey')
ax.axvspan(date(int(2012), 9, 23) ,date(int(2012), 10, 9) , alpha=0.3, color='grey')


ax = axes[1,0]
ax.set_ylabel(r'Mean $V_z$ (mm/d)', fontsize=25)
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
ax.xaxis.set_major_locator(MonthLocator(interval=3))
ax.xaxis.set_minor_locator(MonthLocator(interval=1))
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
fig.autofmt_xdate()
###plot one horizontal line for Vz=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)
###Shade periods corresponding to pumping (26/08/2010 to 15/10/10; 28/09/11 to14/10/11; 23/09/12 to 09/10/12)
ax.axvspan(date(int(2010), 8, 26) ,date(int(2010), 10, 15) , alpha=0.3, color='grey')
ax.axvspan(date(int(2011), 9, 28) ,date(int(2011), 10, 14) , alpha=0.3, color='grey')
ax.axvspan(date(int(2012), 9, 23) ,date(int(2012), 10, 9) , alpha=0.3, color='grey')

ax = axes[1,1]
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
ax.xaxis.set_major_locator(MonthLocator(interval=3))
ax.xaxis.set_minor_locator(MonthLocator(interval=1))
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
fig.autofmt_xdate()
ax.set_title('No Damage', fontsize=26, weight='bold')
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)
###plot one horizontal line for Vz=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)
###Shade periods corresponding to pumping (26/08/2010 to 15/10/10; 28/09/11 to14/10/11; 23/09/12 to 09/10/12)
ax.axvspan(date(int(2010), 8, 26) ,date(int(2010), 10, 15) , alpha=0.3, color='grey')
ax.axvspan(date(int(2011), 9, 28) ,date(int(2011), 10, 14) , alpha=0.3, color='grey')
ax.axvspan(date(int(2012), 9, 23) ,date(int(2012), 10, 9) , alpha=0.3, color='grey')
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
    Step_List = ['Pump2010', 'Refill20102011', 'Pump2011', 'Refill20112012', 'Pump2012', 'Refill20122013']
    StepTsp_List =[1, 5, 1, 5, 1, 5] ##Timestep size (in days) corresponding to simulation step (Step_List)
    StepOut_List =[5, 30, 5, 30, 5, 30] ##Output interval (in days) corresponding to simulation step (Step_List)
    ### Where pressure applies: cavity only, restricted to a conduit, everywhere above cold/temperate transition
    Cases = ['PCavityOnly', 'PRestricted', 'PNotRestric']
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'SigmaI', 'Damage', 'Chi']
    Col_Names_NoD = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'SigmaI']  ##For the no damage (ref) case
    ###START LOOP over each pressure scenario (one subplot par scenario)
    for j, case in enumerate(Cases):
        ### Get the corresponding subplot
        if j ==0:
            ax = axes[0,0]
        elif j ==1:
            ax = axes[0,1]
        elif j==2:
            ax = axes[1,0]
        ax.set_title('{}'.format(case), fontsize=26, weight='bold')
        ###For each pressure scenario get the no damage (ref) case
        ###We create a single dataframe for all steps
        Data_Simu_NoD = pd.DataFrame() ##Initialize an empty dataframe
        if case == 'PRestricted':
            Data_Simu_NoD_WRONG = pd.DataFrame() ##Same for the PRestricted scenario with wrong application of pressure (done as such in Gag sif)
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
            ####For the case PRestricted only, load the version corresponding to wrong application of pressure for comparison (done as such in Gag sif)
            if case =='PRestricted':
                filename_WRONG = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_dm315todm255_Out{}d_WRONG_.dat'.format(case, Step,str(StepTsp),str(StepOut))
                print('Opening file:', filename_WRONG)
                Data_Simu_NoD_WRONG_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_WRONG), names=Col_Names_NoD, delim_whitespace=True)
                ###Drop duplicate lines (bug of SaveLine solver)
                Data_Simu_NoD_WRONG_tmp.drop_duplicates(inplace=True)
                ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
                if i == 0:
                    Data_Simu_NoD_WRONG_tmp['DayOfSimu'] = Data_Simu_NoD_WRONG_tmp['DayOfSimu'] * StepTsp
                else:
                    Data_Simu_NoD_WRONG_tmp['DayOfSimu'] = np.max(Data_Simu_NoD_WRONG['DayOfSimu']) + Data_Simu_NoD_WRONG_tmp['DayOfSimu'] * StepTsp
                data_WRONG = [Data_Simu_NoD_WRONG, Data_Simu_NoD_WRONG_tmp]
                Data_Simu_NoD_WRONG = pd.concat(data_WRONG, ignore_index=True)
        ###Add a column of boolean to dataset depending on wether node is above cavity or not
        Data_Simu_NoD['IsAboveCavity']=IsAboveCavity(Data_Simu_NoD['X'],Data_Simu_NoD['Y'])
        if case == 'PRestricted':
            Data_Simu_NoD_WRONG['IsAboveCavity']=IsAboveCavity(Data_Simu_NoD_WRONG['X'],Data_Simu_NoD_WRONG['Y'])
        ###Store Node above cavity in a separated dataframe
        Data_Simu_NoD_AboveCavity=Data_Simu_NoD[Data_Simu_NoD['IsAboveCavity']]
        if case == 'PRestricted':
            Data_Simu_NoD_WRONG_AboveCavity=Data_Simu_NoD_WRONG[Data_Simu_NoD_WRONG['IsAboveCavity']]
        ###For each timestep calculate mean, min, max of vertical velocity above cavity
        Date = []
        MeanW_NoD = []
        MinW_NoD = []
        MaxW_NoD = []
        if case == 'PRestricted':
            MeanW_NoD_WRONG = []
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
            ###Do the same for the wrong version of PRestricted scenario (mean only)
            if case == 'PRestricted':
                MeanW_NoD_WRONG.append(np.mean(Data_Simu_NoD_WRONG_AboveCavity[Data_Simu_NoD_WRONG_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
        ###Plot MeanW for the case with no damage (ref case) on corresponding subplot
        ax.plot_date(Date, MeanW_NoD, color= 'r', linestyle = '-', linewidth=2, marker='None', xdate=True)
        ###Plot MeanW for the case with no damage (ref case) on subplot 4 for comparison of the 3 pressure scenatio without damage
        if case == 'PCavityOnly':
            Col_PressureScenario='b'
        elif case == 'PRestricted':
            Col_PressureScenario='m'
        elif case == 'PNotRestric':
            Col_PressureScenario='r'
        axes[1,1].plot_date(Date, MeanW_NoD, color=Col_PressureScenario, linestyle='-', linewidth=2, marker='None', xdate=True)
        if case == 'PRestricted':
            axes[1,1].plot_date(Date[0:len(MeanW_NoD_WRONG)], MeanW_NoD_WRONG, color='k', linestyle='--', linewidth=2, marker='None', xdate=True)

        ###NOW DO THE LOOP ON EACH DAMAGE PARAMETER SET
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
                Color = Color_List[l]
                # Color = ColorBlindFriendly_List[l]
            elif Data_ParamSet['Lambdah'][l] == Lambdah_list[1]:
                LineStyle = LineStyle_List[1]
                Color = Color_List[l-8]
                # Color = ColorBlindFriendly_List[l - 8]
            elif Data_ParamSet['Lambdah'][l] == Lambdah_list[2]:
                LineStyle = LineStyle_List[2]
                Color = Color_List[l-16]
                # Color = ColorBlindFriendly_List[l - 16]
            else:
                print('ERROR: Lambdah_list does not correspond to damage parameter sets input file')
            ax.plot_date(Date, MeanW_D, color= Color, linestyle = LineStyle, linewidth=1, marker='None', xdate=True)


    #########################################
    ####   NOW OBSERVATIONS AT STAKES    ####
    #########################################
    ####Load data corresponding to observations at stake
    Pathroot_ObsData = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/PostProcessing/Raw_Data_VeloAtStakes/.')

    ###List of all files containing observationnal data
    ObsDataName_List=['dep14092010-19092010', 'dep19092010-20092010', 'dep20092010-21092010', 'dep21092010-22092010', 'dep22092010-23092010',
    'dep23092010-29092010', 'dep29092010-06102010', 'dep29092010-09092011', 'dep15062011-01072011', 'dep01072011-15072011',
    'dep15072011-10082011', 'dep10082011-31082011', 'dep31082011-09092011', 'dep09092011-23092011', 'dep23092011-28092011',
    'dep28092011-06102011', 'dep06102011-21102011', 'dep21102011-17112011', 'dep06102011-31072012', 'dep13072012-25072012',
    'dep25072012-09082012', 'dep09082012-22082012', 'dep22082012-05092012', 'dep05092012-20092012', 'dep20092012-28092012',
    'dep28092012-01102012', 'dep01102012-03102012']
    ###Name of file containing info on obs data files
    InfoFile_ObsData_Name = 'info_topo2010-2012.dat'
    ###For each year, list of stakes that are not above of the cavity (i.e. out of the AboveCavity mask)
    StakeNotAboveCavity_2010 = [1, 11, 12, 17, 18, 22, 27]
    StakeNotAboveCavity_2011 = [4, 5, 6, 15, 16, 21, 24, 25, 29, 30]
    StakeNotAboveCavity_2012 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 17, 23, 24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
    ###NAMES COLUMS of obs data (same order as in dat.names file)
    Col_Names_ObsData = ['Stake', 'Ux', 'Uy', 'Uz'] ##obs data file contains stake number, displaments x, y and z (m for 2010/2011 and mm for 2012)
    ##The info file contains for each obs data file the two days d0 and d1 between which measurement is performed, the day in between these two days, and the number of days separating these two days
    Col_Names_InfoFile = ['ObsDataFileNumber', 'd0', 'd1', 'dmean', 'NumberOfDays']
    ###LOAD INFO FILE
    InfoFile_ObsData = pd.read_csv(Pathroot_ObsData.joinpath(InfoFile_ObsData_Name), names=Col_Names_InfoFile, delim_whitespace=True)

    ###NOW DO THE LOOP ON ALL OBS DATA FILES
    ###We create a single dataframe for all obs data files
    ObsData = pd.DataFrame() ##Initialize an empty dataframe
    for i,ObsDataName in enumerate(ObsDataName_List):
        year  = ObsDataName[-4:]
        Path_To_ObsFile = 'topo{}/{}.dat'.format(year,ObsDataName)
        ObsData_tmp = pd.read_csv(Pathroot_ObsData.joinpath(Path_To_ObsFile), names=Col_Names_ObsData, delim_whitespace=True)
        ###create additionnal column to dataframe storing, for each stake, info contained in info file
        for k,info in enumerate(Col_Names_InfoFile):
            ObsData_tmp[info]=InfoFile_ObsData[info][i]
        data_obs=[ObsData, ObsData_tmp]
        ObsData=pd.concat(data_obs, ignore_index=True)



    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'RelDensProfile_Forage{}_{}_DInitHL_'.format(ForageNumber,  Case)
    # name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_FullSimu_TestNbIntTsp/'
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()