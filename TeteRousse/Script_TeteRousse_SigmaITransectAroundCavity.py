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

##Function below is used to sort point of contour clockwise
def sort_clockwise(list_of_xy_coords):
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
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)


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
    List_Transect=["AA'","BB'", "CC'"]
    List_Coord_pt1=[[947954.6,2105001.5], [947976.1, 2105120.6], [947943.3, 2105052.1]]
    List_Coord_pt2=[[948027.9,2105114.9], [948024.0,2104971.5], [948112.8,2105049.4]]
    ### Step in the full process from 2010 to 2014
    Step_List = ['P2010', 'Refill20102011']  # , 'P2011', 'Refill20112012', 'P2012', 'Refill20122013']
    StepTsp_List = [1, 5]  # , 1, 5, 1, 5] ##Timestep size (in days) corresponding to simulation step (Step_List)
    StepOut_List = [5, 30]  # , 5, 30, 5, 30] ##Output interval (in days) corresponding to simulation step (Step_List)
    ## Where pressure applies: cavity only: 'PCavityOnly', restricted to a conduit: 'PRestricted', everywhere above cold/temperate transition: 'PNotRestric'
    Case = 'PCav'
    Col_Crevasses_Circ = '#FF0000' ###'#800080' ###Color for representation of circular crevasses
    Col_Crevasses_Other = '#800080' ###Color for representation of non circular crevasses
    Col_Cavity = '#00FFFF'###color for representation of cavity
    Col_Transect = 'lightgrey'###color for representation of transect
    Colormap_for_stressmap = 'cividis' ###Colorm map to choose for map of stress
    ################################################################
    ####  OPEN DATA CORRESPONDING TO CREVASSES AND GROUNDEDMASK ####
    ################################################################
    Pathroot_Obs = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_Crevasses = 'crevasses_xyz_numbered.dat'
    ###load file as dataframe
    Df_Crevasses = pd.read_csv(Pathroot_Obs.joinpath(Filename_Crevasses), names=['X', 'Y', 'Z', 'Crevasse Number', 'IsCircular'], delim_whitespace=True)
    ###Also load glacier contour
    Filename_GlacierContour = 'Contour_TeteRousse2012.dat'
    Df_GlacierContour = pd.read_csv(Pathroot_Obs.joinpath(Filename_GlacierContour), names=['X', 'Y'], delim_whitespace=True)
    xc = Df_GlacierContour['X'].values
    xc=np.append(xc,xc[0]) ##step required to close the contour
    yc = Df_GlacierContour['Y'].values
    yc = np.append(yc, yc[0])  ##step required to close the contour
    ###Open Data corresponding to grounded mask
    Pathroot_GM = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Filename_GM = 'GroundedLine_CLEANED.dat'
    Col_Names_GM = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'GM']
    ###load file as dataframe
    Df_GL = pd.read_csv(Pathroot_GM.joinpath(Filename_GM), names=Col_Names_GM, delim_whitespace=True)
    ###Sort coordinates of cavity contour clockwise
    xy_cavity = np.array([Df_GL['X'].values, Df_GL['Y'].values])
    xy_cavity = np.transpose(xy_cavity)
    xy_cavity_sorted = sort_clockwise(xy_cavity)
    ##close contour
    cavity_contour = np.vstack([xy_cavity_sorted, xy_cavity_sorted[0]])
    #########################################
    ####  OPEN OUTPUT OF THE SIMULATIONS:####
    #########################################
    ### Load files containing surface output of the simulations:
    Pathroot_SimuOutput = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    ###BE WARE: I use the convention SigmaI>SigmaII>SigmaIII which is not the convention of the EigenStress solver in Elmer
    Col_Names_NoD = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Pressure', 'SigmaIII','SigmaII','SigmaI', 'temperature', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz']
    # Col_Names_NoD = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W','SigmaI', 'D','Chi']
    ###We create a single dataframe for all considered steps of the simu
    Data_Simu_NoD = pd.DataFrame()  ##Initialize an empty dataframe
    for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
        filename = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_Out{}d_Tmap_.dat'.format(Case, Step, str(StepTsp), str(StepOut))
        # filename = 'SurfaceOutput_Prono_B097_Sth006_Lh00_{}_{}_tsp{}d_Out{}d_CritMPS_.dat'.format(Case, Step, str(StepTsp), str(StepOut))
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
        ### Create new columns corrresponding to equivalent stress for von Mises, Hayurst and Coulomb criterion
        ##Von Mises:
        Data_Simu_NoD['SigmaVM'] = np.sqrt(0.5 * ((Data_Simu_NoD['SigmaI'] - Data_Simu_NoD['SigmaII']) ** 2 + (Data_Simu_NoD['SigmaI'] - Data_Simu_NoD['SigmaIII']) ** 2 + (Data_Simu_NoD['SigmaII'] - Data_Simu_NoD['SigmaIII']) ** 2))
        ##Hayurst:
        alphaH = 0.21
        betaH = 0.63
        Data_Simu_NoD['SigmaHa'] = alphaH*Data_Simu_NoD['SigmaI']+betaH*Data_Simu_NoD['SigmaVM']+(1-alphaH-betaH)*(-3*Data_Simu_NoD['Pressure'])
        ##Coulomb:
        mu = 0.1
        Data_Simu_NoD['SigmaC'] = 0.5 * (Data_Simu_NoD['SigmaI']-Data_Simu_NoD['SigmaIII']) - 0.5*mu*(Data_Simu_NoD['SigmaI']+Data_Simu_NoD['SigmaIII'])


    # ###################################################################################################
    # ###     PLOT MAP OF EQUIVALENT STRESS AT THE SURFACE AT THE END OF PUMPING FOR EACH CRITERION   ###
    # ###################################################################################################
    ### Get the last day of the step pumping 2010
    Df_LastDayPump2010=Data_Simu_NoD[Data_Simu_NoD['DayOfSimu']==np.max(Data_Simu_NoD[Data_Simu_NoD['Step']=='P2010']['DayOfSimu'])].copy()
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Start loop on the various criterion considered
    for Criterion in ['SigmaI', 'SigmaVM', 'SigmaHa', 'SigmaC']:
        ###Fig 1 is the map of Teterousse, pimp style
        fig1 = plt.figure(1, figsize=(40, 20))
        plt.xlabel(r'X [km]', fontsize=21)
        plt.ylabel(r'Y [km]', fontsize=21)
        plt.tick_params(labelsize=20)  # fontsize of the tick labels
        plt.grid(True)
        plt.grid(alpha=0.5)
        ### Get the figure
        fig1 = plt.figure(1)
        plt.xlim([(np.floor(np.min(Df_LastDayPump2010['X']))-40)/1000, (np.floor(np.max(Df_LastDayPump2010['X']))+40)/1000])
        plt.ylim([(np.floor(np.min(Df_LastDayPump2010['Y']))-40)/1000, (np.floor(np.max(Df_LastDayPump2010['Y']))+40)/1000])
        ###Force x and y axis to same scale
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        ###Create refined regular grid and interpolate field over grid
        x = np.arange(np.floor(np.min(Df_LastDayPump2010['X']))-40,np.floor(np.max(Df_LastDayPump2010['X']))+40,0.1)
        y = np.arange(np.floor(np.min(Df_LastDayPump2010['Y']))-40,np.floor(np.max(Df_LastDayPump2010['Y']))+40,0.1)
        X, Y = np.meshgrid(x, y)
        SigmaEq = Interpolate_field(Df_LastDayPump2010,Criterion,X,Y)
        #shading
        clevs=np.arange(0.0,151,10) ## cbar for shading
        cmap=Colormap_for_stressmap
        #colorbar
        levs_ticks=np.arange(0.0,151,20)
        ###Fills up the map with colors for SigmaEq
        CS1 = plt.contourf(X/1000, Y/1000, SigmaEq*1000, clevs, cmap=cmap,extend='both')
        plt.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ###Show colorbar
        cbar = plt.colorbar(CS1, ticks=levs_ticks, orientation='vertical', label=r'$\sigma_I$ [kPa]')
        ##contour at SigmaEq=0.12 MPa(proposed threshold)
        # CS = plt.contour(X, Y, SigmaEq,[0.0], colors='yellow') ## lines (f(x) contour)
        # plt.clabel(CS, inline=True, fontsize=10) ## numbers in contours
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        ax = plt.gca()
        ax.add_patch(patch)
        for c in CS1.collections:
            c.set_clip_path(patch)
        ####Plot transects over which SigmaEq will be plot on map
        for i,(Transect_Name, Coord_pt1, Coord_pt2) in enumerate(zip(List_Transect, List_Coord_pt1,List_Coord_pt2)):
            ###Use functions to return coords of transect
            coord_transect = Coord_transect(Coord_pt1, Coord_pt2)
            ###Plot transect on map
            plt.plot(coord_transect[0]/1000,coord_transect[1]/1000,color=Col_Transect,linestyle='-',linewidth=4)
            ###Annotate transect name
            if Transect_Name == "AA'":
                plt.annotate(Transect_Name[0], ((coord_transect[0][0]-10)/1000, (coord_transect[1][0]-8)/1000), size=17, weight='bold',color=Col_Transect)
                plt.annotate(Transect_Name[1:], ((coord_transect[0][-1]+3)/1000, (coord_transect[1][-1]-4)/1000), size=17, weight='bold',color=Col_Transect)
            elif Transect_Name == "BB'":
                plt.annotate(Transect_Name[0], ((coord_transect[0][0]-12)/1000, (coord_transect[1][0]+2)/1000), size=17, weight='bold',color=Col_Transect)
                plt.annotate(Transect_Name[1:], ((coord_transect[0][-1]+3)/1000, (coord_transect[1][-1]-9)/1000), size=17, weight='bold',color=Col_Transect)
            elif Transect_Name == "CC'":
                plt.annotate(Transect_Name[0], ((coord_transect[0][0]-12)/1000, (coord_transect[1][0]+2)/1000), size=17, weight='bold',color=Col_Transect)
                plt.annotate(Transect_Name[1:], ((coord_transect[0][-1]+4)/1000, (coord_transect[1][-1]+2)/1000), size=17, weight='bold',color=Col_Transect)
        ###Plot cavity contour
        plt.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.3)
        ###Plot crevasses on map
        # plt.scatter(Df_Crevasses['X']/1000,Df_Crevasses['Y']/1000,color=Col_Crevasses,marker='P',linewidths=0.2)
        ###plot crevasses as continuous line
        for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
            ###get proper points
            Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
            if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
                col = Col_Crevasses_Circ
            else:
                col = Col_Crevasses_Other
            plt.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=1.5)
        ###Show map
        plt.show()
        ###Save map
        # name_output_fig = 'TeteRousseMap_EndPumping2010_Criterion{}_'.format(Criterion)
        # name_path = '/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/PostProcessing/.'
        # path_output_fig = Path(name_path)
        # fig1.savefig(path_output_fig.joinpath(name_output_fig))

    ###########################################################################
    ###     PLOT SIGMA EQUIVALENT OVER CHOSEN TRANSECTS FOR VARIOUS TIME    ###
    ###########################################################################
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Start loop on the various criterion considered
    for Criterion in ['SigmaI', 'SigmaVM', 'SigmaHa', 'SigmaC']:
        ###Fig 2 is subplots for SigmaI profile over transects (row: step of simu, column: transect)
        fig2 = plt.figure(figsize=(60, 60))
        outer = gridspec.GridSpec(2, 3, wspace=0.3, hspace=0.25) ##3 main subplots, 1 on top row, 2 on second row

        # pimp style: transect AA, pumping2010
        ax00 = fig2.add_subplot(outer[0, 0])
        ax00.tick_params(labelsize=18)  # fontsize of the tick labels
        # ax.set_ylim([-9.0, 2.0])
        ax00.grid(True)
        ax00.set_ylabel(r'$\sigma_I$ [kPa]', fontsize=21)
        ax00.yaxis.set_major_formatter(formatter)
        ###plot one horizontal line for Stress=0
        ax00.axhline(y=0.0, color='k', linestyle =':', linewidth=1)

        # pimp style: transect BB, pumping2010:
        ax01 = fig2.add_subplot(outer[0, 1], sharey=ax00)
        # ax.set_ylim([-9.0, 2.0])
        ax01.tick_params(labelsize=18)  # fontsize of the tick labels
        ax01.grid(True)
        ax01.yaxis.set_major_formatter(formatter)
        ###plot one horizontal line for Stress=0
        ax01.axhline(y=0.0, color='k', linestyle =':', linewidth=1)

        # pimp style: transect CC, pumping2010:
        ax02 = fig2.add_subplot(outer[0, 2], sharey=ax00)
        # ax.set_ylim([-9.0, 2.0])
        ax02.tick_params(labelsize=18)  # fontsize of the tick labels
        ax02.grid(True)
        ax02.yaxis.set_major_formatter(formatter)
        ###plot one horizontal line for Stress=0
        ax02.axhline(y=0.0, color='k', linestyle =':', linewidth=1)

        # pimp style: transect AA, Refill2010-2011:
        ax10 = fig2.add_subplot(outer[1, 0], sharex=ax00)
        ax10.tick_params(labelsize=18)  # fontsize of the tick labels
        # ax.set_ylim([-9.0, 2.0])
        ax10.grid(True)
        ax10.set_xlabel(r'Distance [m]', fontsize=21)
        ax10.set_ylabel(r'$\sigma_I$ [kPa]', fontsize=21)
        ax10.yaxis.set_major_formatter(formatter)
        ###plot one horizontal line for Stress=0
        ax10.axhline(y=0.0, color='k', linestyle =':', linewidth=1)

        # pimp style: transect BB, Refill2010-2011:
        ax11 = fig2.add_subplot(outer[1, 1], sharex=ax01, sharey=ax10)
        ax11.set_xlabel(r'Distance [m]', fontsize=21)
        # ax.set_ylim([-9.0, 2.0])
        ax11.tick_params(labelsize=18)  # fontsize of the tick labels
        ax11.grid(True)
        ax11.yaxis.set_major_formatter(formatter)
        ##plot one horizontal line for Stress=0
        ax11.axhline(y=0.0, color='k', linestyle =':', linewidth=1)

        # pimp style: transect CC, Refill2010-2011:
        ax12 = fig2.add_subplot(outer[1, 2], sharex=ax02, sharey=ax10)
        ax12.set_xlabel(r'Distance [m]', fontsize=21)
        # ax.set_ylim([-9.0, 2.0])
        ax12.tick_params(labelsize=18)  # fontsize of the tick labels
        ax12.grid(True)
        ax12.yaxis.set_major_formatter(formatter)
        ##plot one horizontal line for Stress=0
        ax12.axhline(y=0.0, color='k', linestyle =':', linewidth=1)

        #### Start loop over considered transects
        for i,(Transect_Name, Coord_pt1, Coord_pt2) in enumerate(zip(List_Transect, List_Coord_pt1,List_Coord_pt2)):
            ########################################################
            ####  PROCESS CREVASSE DATA TO GET THEM ON TRANSECT ####
            ########################################################
            ###Calculate distance of each crevasse point to transect
            Df_Crevasses['Distance_To_Transect']=Distance_to_line(Coord_pt1, Coord_pt2, [Df_Crevasses['X'], Df_Crevasses['Y']])
            ###Identify the 10 closest crevasses points to the considered transect
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
                        dist = np.sqrt((Df_Crevasses.loc[idx]['X']-Df_Crevasses.loc[indx_ClosestCrevasses[j]]['X'])**2 + (Df_Crevasses.loc[idx]['Y']-Df_Crevasses.loc[indx_ClosestCrevasses[j]]['Y'])**2)
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
            Coord_Crevasses_Projected = Proj_on_line(Coord_pt1,Coord_pt2, [Df_Crevasses.loc[indx_ClosestCrevasses]['X'], Df_Crevasses.loc[indx_ClosestCrevasses]['Y']])
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
                        dist = np.sqrt((Df_GL.loc[idx]['X']-Df_GL.loc[indx_ClosestGL[j]]['X'])**2 + (Df_GL.loc[idx]['Y']-Df_GL.loc[indx_ClosestGL[j]]['Y'])**2)
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
            Coord_GL_Projected = Proj_on_line(Coord_pt1,Coord_pt2, [Df_GL.loc[indx_ClosestGL]['X'], Df_GL.loc[indx_ClosestGL]['Y']])
            ###Convert these coord in terms of distance along transect
            DistOfGL_along_transect=Distance_along_transect(Coord_pt1, Coord_pt2, Coord_GL_Projected)

            # ##ANOTHER METHOD TO COMPUTE CAVITY MASK ALONG TRANSECT
            # ###Interpolate GM on considered transect
            # Interpolated_GM = Interpolate_field(Df_GM,'GM', coord_transect[0], coord_transect[1])
            # ###Find indexes of first non-grounded and last non-grounded nodes
            # for i in range(len(Interpolated_GM)):
            #     if Interpolated_GM[i]>-0.999:
            #         continue
            #     else:
            #         FirstFloating_Idx=i
            #         break
            # for j in range(FirstFloating_Idx,len(Interpolated_GM)):
            #     if Interpolated_GM[j]<-0.999:
            #         continue
            #     else:
            #         LastFloating_Idx=j-1
            #         break
            # DistAlongTransect_CavityStart=dist_along_transect[FirstFloating_Idx]
            # DistAlongTransect_CavityEnd=dist_along_transect[LastFloating_Idx]

            #################################################
            ###NOW GET SIMU RESULTS FOR DAYS OF INTEREST ####
            #################################################
            for Step in Step_List:
                ### We want to product a list corresponding to all simulation days at which we have an output for the considered step
                SimuDays=Data_Simu_NoD[Data_Simu_NoD['Step']==Step]['DayOfSimu'].drop_duplicates()
                if Step=='P2010':###The step corresponding to pumping
                    Ouput_Interval = 5 ##one plot every 5 days
                    ###first get proper subplot
                    if i==0:
                        ax=ax00
                    elif i==1:
                        ax=ax01
                    elif i==2:
                        ax=ax02
                    ax.set_title('Transect {}'.format(Transect_Name), fontsize=21, weight='bold')
                    ###Prepare the ticks for the colorbar
                elif Step=='Refill20102011':
                    Ouput_Interval = 20 ##one plot every 20 days
                    if i==0:
                        ax=ax10
                    elif i==1:
                        ax=ax11
                    elif i==2:
                        ax=ax12
                ###Create a color map to cover all days at which the SigmaI profiles are plot
                cm = plt.get_cmap('RdBu_r')
                cNorm = colors.Normalize(vmin=np.min(SimuDays), vmax=np.max(SimuDays))
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
                ###plot SigmaEq of selected days along considered transect
                k = 0
                levs_ticks_days = [] ##ticks for colorbars
                for j, day in enumerate(SimuDays):
                    ### for pumping step profile every 5 days
                    if not day % Ouput_Interval == 0:
                        continue
                    k = k + 1 ##count total number of days for which a profile is plot
                    levs_ticks_days.append(day)
                    print('Plotting SigmaI profile of step', Step, 'for simu day', day)
                    Data_Simu_NoD_Today = Data_Simu_NoD[Data_Simu_NoD['DayOfSimu'] == day].copy()
                    ###Use functions to return coords of transect
                    coord_transect = Coord_transect(Coord_pt1, Coord_pt2)
                    ###Convert coords of transect in terms of distance along transect
                    dist_along_transect=Distance_along_transect(Coord_pt1, Coord_pt2, coord_transect)
                    ###Interpolate the SigmaEq of the day on the considered transect
                    SigmaEq = Interpolate_field(Data_Simu_NoD_Today, Criterion, coord_transect[0], coord_transect[1])
                    ###Plot interpolated sigmaI as a function of distance along transect
                    ax.plot(dist_along_transect, SigmaEq*1000, color=scalarMap.to_rgba(day), linestyle='-', linewidth=2)
                    ax.set_xlim([-1, np.max(dist_along_transect)+0.5])
                ###plot vertical line corresponding to crevasses
                for l in range(len(DistOfCrevasses_along_transect)):
                    ###a crevasse point can be very close to the line XX' but out of the segment XX' (i.e. out of transect)
                    ###In that case, we do not plot it
                    if DistOfCrevasses_along_transect.values[l] < 0 or DistOfCrevasses_along_transect.values[l]>np.max(dist_along_transect):
                        continue
                    ax.axvline(x=DistOfCrevasses_along_transect.values[l], color='k', linestyle='-', linewidth=3)
                ###Shade area corresponding to cavity based on initial grounded mask
                ax.axvspan(np.min(DistOfGL_along_transect),np.max(DistOfGL_along_transect), alpha=0.3, color='grey')
                if Step == 'P2010' and i == 2:
                    levs_ticks_days.insert(0,1)
                    fig2.colorbar(cmx.ScalarMappable(norm=cNorm, cmap=cm), ax=[ax00, ax01, ax02], ticks=levs_ticks_days[::2], orientation='vertical', label=r'Time [days]')
                elif Step == 'Refill20102011' and i==2:
                    fig2.colorbar(cmx.ScalarMappable(norm=cNorm, cmap=cm), ax=[ax10, ax11, ax12], ticks=levs_ticks_days[::2], orientation='vertical', label=r'Time [days]')
        ### Show Fig:
        plt.show()
        ### Save Fig:
        # name_output_fig2 = 'TeteRousse_ProfilSigmaEqOverTime_3Transects_Criterion{}_'.format(Criterion)
        # name_path = '/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/PostProcessing/.'
        # path_output_fig = Path(name_path)
        # fig2.savefig(path_output_fig.joinpath(name_output_fig2))



