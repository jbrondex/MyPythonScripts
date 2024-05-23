"""
@author: jbrondex

Description:
------------
This file aims at plotting the evolution of maaximum principle stress along some surface transect above the cavity

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
from scipy.interpolate import griddata
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
#####################################################################################
# Define function to interpolate field defined on mesh nodes to the line profile ####
#####################################################################################
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

def Coord_transect(pt1, pt2): ###return list of coordinates corresponding to line passing by pt1 and pt2
    x1, y1 = pt1 ##pt1 defining transect
    x2, y2 = pt2 ##pt2 defining transect
    xmin=min(x1, x2)
    xmax=max(x1, x2)
    x = np.linspace(xmin, xmax, 1000)
    y = ((y2-y1)/(x2-x1))*(x-x1)+y1
    coord = [x,y]
    return coord

def Distance_along_transect(pt1, pt2, pt3): ###return distance of pt3 along transect (pt1,pt2)
    x1, y1 = pt1 ##pt1 defining transect
    x2, y2 =pt2 ##pt2 defining transect
    x, y = pt3
    x0=min(x1, x2)
    if x0 == x1:
        y0 = y1
    elif x0 == x2:
        y0 = y2
    dist_along_line = np.sqrt((x-x0)**2+(y-y0)**2)
    return dist_along_line

def Distance_to_line(pt1, pt2, pt3): ###return distance of pt3 to line define by points (pt1, pt2)
    x1, y1 = pt1 ##pt1 defining transect
    x2, y2 = pt2 ##pt2 defining transect
    x3, y3 = pt3 ##pt3 out of transect. How far is it from transect ?

    det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1) ##This is det(AB;AC) with A and B the two extremities of transect and C is pt3
    distance = det/np.sqrt((x2-x1)**2 + (y2-y1)**2) ##det(AB;AC)/||AB|| = ||AC||*sin(AB,AC) = ||CD||
    return distance

def Proj_on_line(pt1, pt2, pt3):##Return coords of pt3 projected on line defined by (pt1,pt2)
    x1, y1 = pt1
    x2, y2 = pt2
    x3, y3 = pt3
    dx, dy = x2-x1, y2-y1
    scalarprod = dy*(y3-y1)+dx*(x3-x1) ##AB.AC=||AB||*||AC||*cos(AB,AC)=||AB||*||AD||
    a = scalarprod/(dx*dx + dy*dy) ## a = ||AB||*||AD|| / ||AB||**2  = ||AD|| / ||AB||
    return x1+a*dx, y1+a*dy

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
ax.set_ylabel(r'$\sigma_I$ (MPa)', fontsize=25)
ax.tick_params(labelsize=18)  # fontsize of the tick labels
#ax.set_ylim([-9.0, 2.0])
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
###plot one horizontal line for Stress=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)

ax = axes[0,1]
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.set_ylim([-9.0, 2.0])
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
###plot one horizontal line for Stress=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)

ax = axes[1,0]
ax.set_ylabel(r'$\sigma_I$ (MPa)', fontsize=25)
# ax.set_ylim([-9.0, 2.0])
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
###plot one horizontal line for Stress=0
ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)

ax = axes[1,1]
# ax.set_ylim([-9.0, 2.0])
ax.tick_params(labelsize=18)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
ax.set_title('Comp Pressure Scenario', fontsize=26, weight='bold')

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ### General parameter of simu
    RefStartDate = date(int(2011), 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    StartDayNumber_Pumping2010 = "-315"
    StartYear_Pumping2010 = "2010"
    StartDate_Pumping2010 = RefStartDate + timedelta(days=int(StartDayNumber_Pumping2010))
    ###Provide coordinate of points defining lines corresponding to transect on which we want to plot SigmaI: one couple of point per transect
    List_Coord_pt1=[[947954.6,2105001.5], [947976.1, 2105120.6]]
    List_Coord_pt2=[[948027.9,2105114.9], [948024.0,2104971.5]]
    ### Step in the full process from 2010 to 2014
    Step_List = ['Pump2010', 'Refill20102011']  # , 'Pump2011', 'Refill20112012', 'Pump2012', 'Refill20122013']
    StepTsp_List = [1, 5]  # , 1, 5, 1, 5] ##Timestep size (in days) corresponding to simulation step (Step_List)
    StepOut_List = [5, 30]  # , 5, 30, 5, 30] ##Output interval (in days) corresponding to simulation step (Step_List)
    ## Where pressure applies: cavity only: 'PCavityOnly', restricted to a conduit: 'PRestricted', everywhere above cold/temperate transition: 'PNotRestric'
    Case = 'PRestricted'

    ################################################################
    ####  OPEN DATA CORRESPONDING TO CREVASSES AND GROUNDEDMASK ####
    ################################################################
    Pathroot_Obs = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_Crevasses = 'crevasses_xyz.dat'
    ###load file as dataframe
    Df_Crevasses = pd.read_csv(Pathroot_Obs.joinpath(Filename_Crevasses), names=['X', 'Y', 'Z'], delim_whitespace=True)
    ###Open Data corresponding to grounded mask
    Pathroot_GM = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Filename_GM = 'GroundedMaskInit.dat'
    Col_Names_GM = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'GM']
    ###load file as dataframe
    Df_GM = pd.read_csv(Pathroot_GM.joinpath(Filename_GM), names=Col_Names_GM, delim_whitespace=True)
    Df_GM.drop_duplicates(inplace=True)
    ###Keep only the Grounding Line
    Df_GL=Df_GM[Df_GM['GM']==0]

    #########################################
    ####  OPEN OUTPUT OF THE SIMULATIONS:####
    #########################################
    ### Load files containing surface output of the simulations:
    Pathroot_SimuOutput = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names_NoD = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W','SigmaI']
    ###We create a single dataframe for all considered steps of the simu
    Data_Simu_NoD = pd.DataFrame()  ##Initialize an empty dataframe
    for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
        filename = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_dm315todm255_Out{}d_.dat'.format(Case, Step, str(StepTsp), str(StepOut))
        print('Opening file:', filename)
        Data_Simu_NoD_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names_NoD, delim_whitespace=True)
        ###Drop duplicate lines (bug of SaveLine solver)
        Data_Simu_NoD_tmp.drop_duplicates(inplace=True)
        ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
        if i == 0:
            Data_Simu_NoD_tmp['DayOfSimu'] = Data_Simu_NoD_tmp['DayOfSimu'] * StepTsp
        else:
            Data_Simu_NoD_tmp['DayOfSimu'] = np.max(Data_Simu_NoD['DayOfSimu']) + Data_Simu_NoD_tmp['DayOfSimu'] * StepTsp
        ###We create an additionnal column storing the step of simu we are dealing with
        Data_Simu_NoD_tmp['Step'] = Step
        data = [Data_Simu_NoD, Data_Simu_NoD_tmp]
        Data_Simu_NoD = pd.concat(data, ignore_index=True)
        ### We want to product a list corresponding to all simulation days at which we have an output
    SimuDays = Data_Simu_NoD['DayOfSimu']
    SimuDays.drop_duplicates(inplace=True)
    ###Create a color map to cover all combinations of (Sigmath, B) for a given value of lamddah
    cm = plt.get_cmap('RdBu_r')
    cNorm = colors.Normalize(vmin=0, vmax=41)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    Color_List = []
    for i in range(41):
        Color_List.append(scalarMap.to_rgba(i))

    ##################################################
    ####  DO A LOOP OVER EACH CONSIDERED TRANSECT ####
    ##################################################
    for i,(Coord_pt1, Coord_pt2) in enumerate(zip(List_Coord_pt1,List_Coord_pt2)):
        ### Get the corresponding subplot
        if i == 0:
            ax = axes[0, 0]
        elif i == 1:
            ax = axes[0, 1]
        # elif j == 2:
        #     ax = axes[1, 0]
        ax.set_title('Transect n° {}'.format(str(i+1)), fontsize=26, weight='bold')

        ###Use functions to return coords of transect
        coord_transect = Coord_transect(Coord_pt1, Coord_pt2)
        ###Convert coords of transect in terms of distance along transect
        dist_along_transect=Distance_along_transect(Coord_pt1, Coord_pt2, coord_transect)
        ########################################################
        ####  PROCESS CREVASSE DATA TO GET THEM ON TRANSECT ####
        ########################################################
        ###Calculate distance of each crevasse point to transect
        Df_Crevasses['Distance_To_Transect']=Distance_to_line(Coord_pt1, Coord_pt2, [Df_Crevasses['X'], Df_Crevasses['Y']])
        ###Identify the 10 closest crevasses to the considered transect
        Ten_Closest_Crevasses = Df_Crevasses['Distance_To_Transect'].abs().nsmallest(10)
        ### Among the 10 closest crevasses points identified above we want to remove those belonging to the same crevasse
        ### For that we eliminate points that are too close (i.e. below a distance threshold) of points already identified
        indx_ClosestCrevasses=[]
        for k, idx in enumerate(Ten_Closest_Crevasses.index):
            if k==0:
                indx_ClosestCrevasses.append(idx)
            else:
                for j in range(len(indx_ClosestCrevasses)): ##Look at the distance of current point to all those that have been already kept
                    print('Comparing pt',idx, 'to pt',indx_ClosestCrevasses[j] )
                    dist = np.sqrt((Df_Crevasses.iloc[idx]['X']-Df_Crevasses.iloc[indx_ClosestCrevasses[j]]['X'])**2 + (Df_Crevasses.iloc[idx]['Y']-Df_Crevasses.iloc[indx_ClosestCrevasses[j]]['Y'])**2)
                    print('dist is', dist)
                    if dist < 5: ##Points belong to same crevasse
                        print('we discard pt', idx)
                        break
                    else:
                        if j == len(indx_ClosestCrevasses)-1:
                            print('we keep pt', idx)
                            indx_ClosestCrevasses.append(idx)
                        else:
                            continue
        ###Now project the coordinates of crevasses on the considered transect
        Coord_Crevasses_Projected = Proj_on_line(Coord_pt1,Coord_pt2, [Df_Crevasses.iloc[indx_ClosestCrevasses]['X'], Df_Crevasses.iloc[indx_ClosestCrevasses]['Y']])
        ###Convert these coord in terms of distance along transect
        DistOfCrevasses_along_transect=Distance_along_transect(Coord_pt1, Coord_pt2, Coord_Crevasses_Projected)

        ###############################################
        ###NOW PROCESS DATA REGARDING GROUNDEDMASK ####
        ###############################################
        ###Calculate distance of each point of the GL to transect
        DistFromGLToTransect=Distance_to_line(Coord_pt1, Coord_pt2, [Df_GL['X'], Df_GL['Y']])
        ###Identify the 5 closest GL points to the considered transect
        Five_Closest_GL = DistFromGLToTransect.abs().nsmallest(5)
        ### Among the 5 closest GL points identified above we want to identify the two extremities of the cavity
        ### For that we eliminate GL points that are too close (i.e. below a distance threshold) of points already identified
        indx_ClosestGL=[]
        for k, idx in enumerate(Five_Closest_GL.index):
            if k==0:
                indx_ClosestGL.append(idx)
            else:
                for j in range(len(indx_ClosestGL)): ##Look at the distance of current point to all those that have been already kept
                    print('Comparing pt',idx, 'to pt',indx_ClosestGL[j] )
                    dist = np.sqrt((Df_GM.iloc[idx]['X']-Df_GM.iloc[indx_ClosestGL[j]]['X'])**2 + (Df_GM.iloc[idx]['Y']-Df_GM.iloc[indx_ClosestGL[j]]['Y'])**2)
                    print('dist is', dist)
                    if dist < 40: ##Points belong to same crevasse
                        print('we discard pt', idx)
                        break
                    else:
                        if j == len(indx_ClosestGL)-1:
                            print('we keep pt', idx)
                            indx_ClosestGL.append(idx)
                        else:
                            continue
        ###Now project the coordinates of crevasses on the considered transect
        Coord_GL_Projected = Proj_on_line(Coord_pt1,Coord_pt2, [Df_GM.iloc[indx_ClosestGL]['X'], Df_GM.iloc[indx_ClosestGL]['Y']])
        ###Convert these coord in terms of distance along transect
        DistOfGL_along_transect=Distance_along_transect(Coord_pt1, Coord_pt2, Coord_GL_Projected)

        ##ANOTHER METHOD TO COMPUTE CAVITY MASK ALONG TRANSECT
        ###Interpolate GM on considered transect
        Interpolated_GM = Interpolate_field(Df_GM,'GM', coord_transect[0], coord_transect[1])
        ###Find indexes of first non-grounded and last non-grounded nodes
        for i in range(len(Interpolated_GM)):
            if Interpolated_GM[i]>-0.999:
                continue
            else:
                FirstFloating_Idx=i
                break
        for j in range(FirstFloating_Idx,len(Interpolated_GM)):
            if Interpolated_GM[j]<-0.999:
                continue
            else:
                LastFloating_Idx=j-1
                break
        DistAlongTransect_CavityStart=dist_along_transect[FirstFloating_Idx]
        DistAlongTransect_CavityEnd=dist_along_transect[LastFloating_Idx]

        ###plot SigmaI of selected days along considered transect
        k = 0
        for j, day in enumerate(SimuDays):
            ### for now one profile every 10 days
            if not day % 10 == 0:
                continue
            k = k + 1 ##count total number of days for which a profile is plot
            print('Plotting SigmaI profile for simu day', day)
            Data_Simu_NoD_Today = Data_Simu_NoD[Data_Simu_NoD['DayOfSimu'] == day]
            ###Interpolate the SigmaI of the day on the considered transect
            Interpolated_SigmaI = Interpolate_field(Data_Simu_NoD_Today, 'SigmaI', coord_transect[0], coord_transect[1])
            ###Plot interpolated sigmaI as a function of distance along transect
            ax.plot(dist_along_transect, Interpolated_SigmaI, color=Color_List[k], linestyle='-', linewidth=2)
        ###plot vertical line corresponding to crevasses
        for l in range(len(DistOfCrevasses_along_transect)):
            ax.axvline(x=DistOfCrevasses_along_transect.values[l], color='k', linestyle='-', linewidth=3)
        ###Shade area corresponding to cavity based on initial grounded mask
        ax.axvspan(np.min(DistOfGL_along_transect),np.max(DistOfGL_along_transect), alpha=0.3, color='grey')
        ax.axvspan(DistAlongTransect_CavityStart,DistAlongTransect_CavityEnd, alpha=0.3, color='red')




    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'RelDensProfile_Forage{}_{}_DInitHL_'.format(ForageNumber,  Case)
    # name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_FullSimu_TestNbIntTsp/'
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()