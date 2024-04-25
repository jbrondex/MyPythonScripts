"""
@author: jbrondex

Description:
------------
This file makes subplots at different times of density profiles obtained at a given forage of the Taconnaz glacier with the Particle Advector solver.
The purpose is to evaluate the impact of the number of internal timesteps on the solutions produced by the Particle solver.
Here, simulations are run with true velocity (Porous is solved) and there is no analytical solution.

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

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

fig, axes = plt.subplots(1, 4, figsize=(60, 20))
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

ax = axes[3]
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
    Case = 'Full Simu' #'UniformVelo'
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'Year', 'ForageNumber', 'NodeNumber', 'X', 'Y', 'Z', 'Depth', 'DensRel', 'Porous 3']
    ###List all combinations of Method and TimeStep to plot on a subplot (with corresponding line style/color)
    IntTsp_List=[1, 5, 10, 20]
    LineStyle_List = ['-', '-.', '--', ':']
    #### Color List from a color map
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=3)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    Color_List=[]
    for i in range(np.size(IntTsp_List)):
        Color_List.append(scalarMap.to_rgba(i))

    ### Years of interest (one subplot per year):
    Year_List =[1, 5 , 15, 79]

    ### GO FOR IT:
    ### Give the name of Forage Number to the Figure
    Name_Fig = 'Forage N°{}, Case {}'.format(ForageNumber, Case)
    fig.suptitle(Name_Fig, fontsize=36, weight='bold')
    #### Load data for all combination in a DataFrame
    for l, (IntTsp, LineCol, LineSty) in enumerate(zip(IntTsp_List, Color_List, LineStyle_List)):
        print(l)
        Data=[]
        ### Open file corresponding to each partition and then concatenate in single DataFrame
        for part in range(Npartition):
            if IntTsp == 10:
                file_name='forage_MyTaco_Porous_SL_tsp005_Out1a_80a_FromSteady_NoBC5a8PTgt_NoPWallTop_DensBF_.dat.{}'.format(part)
            else:
                file_name = 'forage_MyTaco_Porous_SL_tsp005_Out1a_80a_FromSteady_NoBC5a8PTgt_NoPWallTop_DensBF_{}IntTsp_.dat.{}'.format(IntTsp, part)
            file = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names, delim_whitespace=True)
            Data.append(file)
        Data=pd.concat(Data)
        #### Keep only data for the forage of interest and in the ice body (DensRel > 0)
        Data_Forage = Data[(Data['ForageNumber']==ForageNumber) & (Data['DensRel']>0)]
        #### Plot current combination for each of the considered years
        for k,year in enumerate(Year_List):
            ax=axes[k]
            ax.set_title('@{}a'.format(year), fontsize=34, weight='bold')
            ### Plot initial D profile (Profile of H and L) -> ONCE ONLY
            if l==0:
                DensityInit = np.maximum(1.0 - 0.55 * np.exp(-0.038 * Data_Forage[Data_Forage['Year'] == 1]['Depth'].values), 350/917)
                DensityInit[-1] = 350/917 ### Top initial condition on RelDens corresponds to BC condition
                ax.plot(DensityInit, - Data_Forage[Data_Forage['Year'] == 1]['Depth'].values, color= 'k', linestyle = '-', linewidth=2)
            year = year + 1
            if Data_Forage[Data_Forage['Year']==year]['DensRel'].empty:  ####If no results for the considered year (simulation still running) then go to next case
                continue
            print('Number of IntTsp :', IntTsp, 'Time Step n° ',Data_Forage[Data_Forage['Year']==year]['Timestep'].values[0], ', LineStyle = ', LineSty)
            ###MESH: Plot horizontal lines corresponding to node position (in solid for DensSol and dashed for SL)
            for i in range(len(Data_Forage[Data_Forage['Year'] == year]['Depth'].values)):
                ax.axhline(y=-Data_Forage[Data_Forage['Year'] == year]['Depth'].values[i], color='darkgray',linestyle='--', linewidth=0.8)
            #### Now plot density profiles:
            ax.plot(Data_Forage[Data_Forage['Year']==year]['DensRel'].values, - Data_Forage[Data_Forage['Year']==year]['Depth'].values, color= LineCol, linestyle = LineSty, linewidth=4.5)
            ax.set_ylim([np.min(- Data_Forage[Data_Forage['Year']==year]['Depth']) - 5, 0])

    #### DUMMY plot for legend
    ax = axes[1]
    for m in range(np.size(IntTsp_List)):
        print(m)
        Label = '{} Int Tsp'.format(IntTsp_List[m])
        print(Label)
        ax.plot(np.NaN, np.NaN, label=Label, color=Color_List[m], linewidth=4, linestyle=LineStyle_List[m])
    fig.subplots_adjust(bottom=0.15, wspace=0.25)
    fig.legend(loc='lower center', fancybox=True, shadow=True,  fontsize=30, ncol=2)


    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'RelDensProfile_Forage{}_{}_DInitHL_'.format(ForageNumber,  Case)
    # name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_FullSimu_TestNbIntTsp/'
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()