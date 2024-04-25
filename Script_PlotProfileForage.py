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

fig, axes = plt.subplots(1, 3, figsize=(40, 20))
# pimp style
##First Row
ax = axes[0]
ax.set_ylabel(r'Depth (m)', fontsize=34)
ax.set_xlabel(r'$\rho / \rho_i$', fontsize=34)
ax.set_xlim([0.3, 1.1])
ax.tick_params(labelsize=24)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
###Plot vertical lines corresponding to min and max rel density
ax.axvline(x=0.38, color= 'k', linestyle = ':', linewidth=3)
ax.axvline(x=1, color= 'k', linestyle = ':', linewidth=3)

ax = axes[1]
ax.set_xlabel(r'$\rho / \rho_i$', fontsize=34)
ax.set_xlim([0.3, 1.1])
ax.tick_params(labelsize=24)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
###Plot vertical lines corresponding to min and max rel density
ax.axvline(x=0.38, color= 'k', linestyle = ':', linewidth=3)
ax.axvline(x=1, color= 'k', linestyle = ':', linewidth=3)

ax = axes[2]
ax.set_xlabel(r'$\rho / \rho_i$', fontsize=34)
ax.set_xlim([0.3, 1.1])
ax.tick_params(labelsize=24)  # fontsize of the tick labels
# ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
###Plot vertical lines corresponding to min and max rel density
ax.axvline(x=0.38, color= 'k', linestyle = ':', linewidth=3)
ax.axvline(x=1, color= 'k', linestyle = ':', linewidth=3)



################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/ForageOutput/.')
    ### Number of Partitions
    Npartition = 4
    ### Number of the forage of interest (FROM 1 TO 11):
    ForageNumber = 11
    ### Considered case : ConstantVelo (i.e. velo field from steady state) or UniformVelo (Porous 1 = Porous 2 =0, Porous 3 = 3 m/a)
    Case = 'UniformVelo' #'UniformVelo'
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'Year', 'ForageNumber', 'NodeNumber', 'X', 'Y', 'Z', 'Depth', 'DensRel', 'Porous 3']
    ###List all combinations of Method and TimeStep to plot on a subplot (with corresponding line style/color)
    Solver_list = ['DensSolv', 'DensSolv', 'DensSolv', 'SL', 'SL', 'SL']
    TspSize_list = ['tsp1', 'tsp01', 'tsp001', 'tsp1', 'tsp01', 'tsp001']
    Color_List = ['limegreen', 'limegreen', 'limegreen', 'darkred', 'darkred', 'darkred']
    LineStyle_List = ['-', '--', ':', '-', '--', ':']
    ### Years of interest (one subplot per year):
    Year_List =[1, 10 ,19]

    ### GO FOR IT:
    ### Give the name of Forage Number to the Figure
    Name_Fig = 'Forage N°{}, Case {}'.format(ForageNumber, Case)
    fig.suptitle(Name_Fig, fontsize=36, weight='bold')
    #### Load data for all combination in a DataFrame
    for (solver, dt, LineCol, LineSty) in zip(Solver_list, TspSize_list, Color_List, LineStyle_List):
        Label = '{}_{}'.format(solver, dt)
        Data=[]
        ### Open file corresponding to each partition and then concatenate in single DataFrame
        for part in range(Npartition):
            if solver == 'DensSolv':
                file_name='forage_MyTaco_Porous_DensSolv_{}_20a_Out1a_{}_NoSource_DInitHL_.dat.{}'.format(dt,  Case, part)
            elif solver == 'SL':
                file_name = 'forage_MyTaco_Porous_SL_{}_20a_Out1a_{}_BC5And8PartTgt_NoSource_DInitHL_.dat.{}'.format(dt, Case, part)
            file = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names, delim_whitespace=True)
            Data.append(file)
        Data=pd.concat(Data)
        #### Keep only data for the forage of interest and in the ice body (DensRel > 0)
        Data_Forage = Data[(Data['ForageNumber']==ForageNumber) & (Data['DensRel']>0)]
        ###MESH: Plot horizontal lines corresponding to node position
        for i in range(len(Data_Forage[Data_Forage['Year'] == 2]['Depth'].values)):
            axes[0].axhline(y=-Data_Forage[Data_Forage['Year'] == 2]['Depth'].values[i], color='darkgray', linestyle='-',linewidth=0.8)
            axes[1].axhline(y=-Data_Forage[Data_Forage['Year'] == 2]['Depth'].values[i], color='darkgray', linestyle='-',linewidth=0.8)
            axes[2].axhline(y=-Data_Forage[Data_Forage['Year'] == 2]['Depth'].values[i], color='darkgray', linestyle='-', linewidth=0.8)
        #### Plot current combination for each of the considered years
        for k,year in enumerate(Year_List):
            ax=axes[k]
            ax.set_title('@{}a'.format(year), fontsize=34, weight='bold')
            if Case == 'UniformVelo': ##In that case we have an analytic solution, show it on graph !
                DepthAnalytique = np.linspace(0,150,1000)
                DensityAnalytique = np.maximum(1.0-0.55*np.exp(-0.038*(DepthAnalytique- 3 * year)), 350/917) #Density(depth) at year = Density(depth - distance parcourue) at initial time
                ax.plot(DensityAnalytique, - DepthAnalytique, color= 'k', linestyle = '-', linewidth=2)
            if not dt == 'tsp1': ##Because first time Save Data solver is called is after one timetep not after one year
                year = year + 1
            if Data_Forage[Data_Forage['Year']==year]['DensRel'].empty:  ####If no results for the considered year (simulation still running) then go to next case
                continue
            print('Solver :', solver,', Time Step size = ', dt,'a', ', Time Step n° ',Data_Forage[Data_Forage['Year']==year]['Timestep'].values[0])
            ax.plot(Data_Forage[Data_Forage['Year']==year]['DensRel'].values, - Data_Forage[Data_Forage['Year']==year]['Depth'].values, color= LineCol, linestyle = LineSty, linewidth=4.5)
            ax.set_ylim([np.min(- Data_Forage[Data_Forage['Year']==year]['Depth']) - 5, 0])

    #### DUMMY plot for legend
    ax = axes[1]
    ax.plot(np.NaN, np.NaN, label='Dens Solv', color='limegreen', linewidth=4, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='Semi Lag', color='darkred', linewidth=4, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =1~\mathrm{a}$', color='darkgrey', linewidth=4, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =0.1~\mathrm{a}$', color='darkgrey', linewidth=4, linestyle='--')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =0.01~\mathrm{a}$', color='darkgrey', linewidth=4, linestyle=':')
    fig.subplots_adjust(bottom=0.15, wspace=0.25)
    fig.legend(loc='lower center', fancybox=True, shadow=True,  fontsize=30, ncol=4)


    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    name_output_fig = 'RelDensProfile_Forage{}_Case{}_DInitHL_'.format(ForageNumber,  Case)
    name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_{}_DInitHL/.'.format(Case)
    path_output_fig = Path(name_path)
    fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()