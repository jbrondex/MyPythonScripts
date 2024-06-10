"""
@author: jbrondex

Description:
------------
This file makes subplots of horizontal and vertical profiles of the three principle stresses obtained with the Stokes solver with linear visco on the one hand,
and the linear elastic solver (called stress solver in Elmer) with various young modulus and poisson coefficient on the other hand. When the poisson coefficient
tends towards 0.5 (incompressible) the stress solutions obtained assuming linear elasticity should tend towards the viscous stresses
"""


################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt

from pathlib import Path

import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files

from matplotlib import ticker
import matplotlib.colors as colors
import matplotlib.cm as cmx

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

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

###Parameters if figure not saved
fig2, axes = plt.subplots(2, 3, figsize=(60, 60))
fig2.subplots_adjust(hspace=0.2,wspace=0.35)
for i in range(0,2):
    for j in range(0,3):
        ax = axes[i,j]
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        ax.grid(True)
        ax.yaxis.set_major_formatter(formatter)
        if i==0: ##Top Line Horizontal profile
            ax.set_xlim([0.0, 10.0])
            ax.set_xlabel(r'Distance [m]', fontsize=21)
            if j==0:
                ax.set_ylabel(r'$\sigma$ [MPa]', fontsize=21)
        elif i==1: ##Bottom Line Vertical profile
            ax.set_ylim([0.0, 10.0])
            ax.set_xlabel(r'$\sigma$ [MPa]', fontsize=21)
            if j==0:
                ax.set_ylabel(r'Distance [m]', fontsize=21)

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/Documents/TestCase_LinearElasticVSStokesLine/ScalarOutput/.')

    ### Load the viscous case (reference !!)
    ##BE AWARE: the eigen stress are not in the same order for elastic and viscous (eigenstress 1 and 3 for the viscous are resp. min and max)
    Col_Names_Viscous = ['Tsp','Line','Node','X','Y','Z','SigmaMin','SigmaMiddle','SigmaMax']
    file_name='Output_Test_StokesLine_BCBottomLeftNode_.dat'
    df_Viscous = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names_Viscous, delim_whitespace=True)
    ###Same for non-linear case
    file_name_NL='Output_Test_StokesNonLine_BCBottomLeftNode_.dat'
    df_Viscous_NL = pd.read_csv(pathroot_mycode.joinpath(file_name_NL), names=Col_Names_Viscous, delim_whitespace=True)
    ### For the elastic case
    ###List all combinations of parameters in the case of the elastic solver
    ##BE AWARE: the eigen stress are not in the same order for elastic and viscous (Principal Stress 1 and 3 for the elastic are resp. max and min)
    Col_Names_Elastic = ['Tsp','Line','Node','X','Y','Z','SigmaMax','SigmaMiddle','SigmaMin']
    PoissonName_List = ['030', '035', '040', '045', '049', '04999']
    Poisson_List =[0.3,0.35,0.4,0.45,0.49,0.4999]
    Young_List = ['1', '9']
    LineStyle_List = [':', '-']

    ### Start the loop on the horizontal (Line 1) and Vertical (Line 2) profiles (= rows of subplot)
    for i,Line in enumerate([1,2]):
        ### We select only outputs along the current line
        df_Viscous_tmp=df_Viscous[df_Viscous['Line']==Line]
        df_Viscous_NL_tmp=df_Viscous_NL[df_Viscous_NL['Line']==Line]
        ## Loop on the 3 principle stresses (columns of subplot)
        for j,Stress in enumerate(['SigmaMax', 'SigmaMiddle', 'SigmaMin']):
            ##Get proper ax
            ax=axes[i,j]
            ##Add title to columns
            if i==0:
                if j==0:
                    ax.set_title(r'$\sigma_I$', fontsize=21, weight='bold')
                elif j==1:
                    ax.set_title(r'$\sigma_{II}$', fontsize=21, weight='bold')
                else:
                    ax.set_title(r'$\sigma_{III}$', fontsize=21, weight='bold')
            ##Loop on combinations of elasticity
            ###Create a color map to cover all poisson number tested
            cm = plt.get_cmap('viridis_r')
            cNorm = colors.Normalize(vmin=0.3, vmax=0.5)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
            for k,(PoissonName,PoissonCoeff) in enumerate(zip(PoissonName_List,Poisson_List)):
                for l,(Young,LineStyle) in enumerate(zip(Young_List, LineStyle_List)):
                    filename = 'Output_Test_Elastic_BCBottomLeftNode_Poisson{}_Young{}GPa_.dat'.format(PoissonName, Young)
                    print('Opening file:', filename)
                    df_Elastic = pd.read_csv(pathroot_mycode.joinpath(filename), names=Col_Names_Elastic, delim_whitespace=True)
                    df_Elastic_tmp=df_Elastic[df_Elastic['Line']==Line]
                    ##Now plot the current elastic stress profile
                    if i==0:
                        df_Elastic_tmp.sort_values(by=['X'], inplace=True)
                        ax.plot(df_Elastic_tmp['X'].values,df_Elastic_tmp[Stress].values, color=scalarMap.to_rgba(PoissonCoeff), linewidth=2, linestyle=LineStyle)
                    else:
                        df_Elastic_tmp.sort_values(by=['Y'], inplace=True)
                        ax.plot(df_Elastic_tmp[Stress].values,df_Elastic_tmp['Y'].values, color=scalarMap.to_rgba(PoissonCoeff), linewidth=2, linestyle=LineStyle)
            ##Plot the reference stress profile (i.e. viscous stress) on top
            if i == 0:
                df_Viscous_tmp.sort_values(by=['X'], inplace=True)
                df_Viscous_NL_tmp.sort_values(by=['X'], inplace=True)
                ax.plot(df_Viscous_tmp['X'].values, df_Viscous_tmp[Stress].values, color='r', linewidth=2.2,linestyle='-')
                ax.plot(df_Viscous_NL_tmp['X'].values, df_Viscous_NL_tmp[Stress].values, color=SkyBlue, linewidth=2.2,linestyle='-')
            else:
                df_Viscous_tmp.sort_values(by=['Y'], inplace=True)
                df_Viscous_NL_tmp.sort_values(by=['Y'], inplace=True)
                ax.plot(df_Viscous_tmp[Stress].values, df_Viscous_tmp['Y'].values, color='r', linewidth=2.2,linestyle='-')
                ax.plot(df_Viscous_NL_tmp[Stress].values, df_Viscous_NL_tmp['Y'].values, color=SkyBlue, linewidth=2.2,linestyle='-')

    fig2.colorbar(cmx.ScalarMappable(norm=cNorm, cmap=cm), ax=axes[0,:], ticks=[0.3,0.35,0.4,0.45,0.5],orientation='vertical', label=r'Poisson Coeff')
    fig2.colorbar(cmx.ScalarMappable(norm=cNorm, cmap=cm), ax=axes[1,:], ticks=[0.3,0.35,0.4,0.45,0.5],orientation='vertical', label=r'Poisson Coeff')
    plt.show()





    # ### GO FOR IT:
    # ### Give the name of Forage Number to the Figure
    # Name_Fig = 'Forage N°{}, Case {}'.format(ForageNumber, Case)
    # fig.suptitle(Name_Fig, fontsize=36, weight='bold')
    # #### Load data for all combination in a DataFrame
    # for (solver, dt, LineCol, LineSty) in zip(Solver_list, TspSize_list, Color_List, LineStyle_List):
    #     Label = '{}_{}'.format(solver, dt)
    #     Data=[]
    #     ### Open file corresponding to each partition and then concatenate in single DataFrame
    #     for part in range(Npartition):
    #         if solver == 'DensSolv':
    #             file_name='forage_MyTaco_Porous_DensSolv_{}_20a_Out1a_{}_NoSource_DInitHL_.dat.{}'.format(dt,  Case, part)
    #         elif solver == 'SL':
    #             file_name = 'forage_MyTaco_Porous_SL_{}_20a_Out1a_{}_BC5And8PartTgt_NoSource_DInitHL_.dat.{}'.format(dt, Case, part)
    #         file = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names, delim_whitespace=True)
    #         Data.append(file)
    #     Data=pd.concat(Data)
    #     #### Keep only data for the forage of interest and in the ice body (DensRel > 0)
    #     Data_Forage = Data[(Data['ForageNumber']==ForageNumber) & (Data['DensRel']>0)]
    #     ###MESH: Plot horizontal lines corresponding to node position
    #     for i in range(len(Data_Forage[Data_Forage['Year'] == 2]['Depth'].values)):
    #         axes[0].axhline(y=-Data_Forage[Data_Forage['Year'] == 2]['Depth'].values[i], color='darkgray', linestyle='-',linewidth=0.8)
    #         axes[1].axhline(y=-Data_Forage[Data_Forage['Year'] == 2]['Depth'].values[i], color='darkgray', linestyle='-',linewidth=0.8)
    #         axes[2].axhline(y=-Data_Forage[Data_Forage['Year'] == 2]['Depth'].values[i], color='darkgray', linestyle='-', linewidth=0.8)
    #     #### Plot current combination for each of the considered years
    #     for k,year in enumerate(Year_List):
    #         ax=axes[k]
    #         ax.set_title('@{}a'.format(year), fontsize=34, weight='bold')
    #         if Case == 'UniformVelo': ##In that case we have an analytic solution, show it on graph !
    #             DepthAnalytique = np.linspace(0,150,1000)
    #             DensityAnalytique = np.maximum(1.0-0.55*np.exp(-0.038*(DepthAnalytique- 3 * year)), 350/917) #Density(depth) at year = Density(depth - distance parcourue) at initial time
    #             ax.plot(DensityAnalytique, - DepthAnalytique, color= 'k', linestyle = '-', linewidth=2)
    #         if not dt == 'tsp1': ##Because first time Save Data solver is called is after one timetep not after one year
    #             year = year + 1
    #         if Data_Forage[Data_Forage['Year']==year]['DensRel'].empty:  ####If no results for the considered year (simulation still running) then go to next case
    #             continue
    #         print('Solver :', solver,', Time Step size = ', dt,'a', ', Time Step n° ',Data_Forage[Data_Forage['Year']==year]['Timestep'].values[0])
    #         ax.plot(Data_Forage[Data_Forage['Year']==year]['DensRel'].values, - Data_Forage[Data_Forage['Year']==year]['Depth'].values, color= LineCol, linestyle = LineSty, linewidth=4.5)
    #         ax.set_ylim([np.min(- Data_Forage[Data_Forage['Year']==year]['Depth']) - 5, 0])
    #
    # #### DUMMY plot for legend
    # ax = axes[1]
    # ax.plot(np.NaN, np.NaN, label='Dens Solv', color='limegreen', linewidth=4, linestyle='-')
    # ax.plot(np.NaN, np.NaN, label='Semi Lag', color='darkred', linewidth=4, linestyle='-')
    # ax.plot(np.NaN, np.NaN, label='$\Delta t =1~\mathrm{a}$', color='darkgrey', linewidth=4, linestyle='-')
    # ax.plot(np.NaN, np.NaN, label='$\Delta t =0.1~\mathrm{a}$', color='darkgrey', linewidth=4, linestyle='--')
    # ax.plot(np.NaN, np.NaN, label='$\Delta t =0.01~\mathrm{a}$', color='darkgrey', linewidth=4, linestyle=':')
    # fig.subplots_adjust(bottom=0.15, wspace=0.25)
    # fig.legend(loc='lower center', fancybox=True, shadow=True,  fontsize=30, ncol=4)
    #
    #
    # ################################################################################
    # # SAVE THE FIGURES #####
    # ################################################################################
    # name_output_fig = 'RelDensProfile_Forage{}_Case{}_DInitHL_'.format(ForageNumber,  Case)
    # name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_{}_DInitHL/.'.format(Case)
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

