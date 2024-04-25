"""
@author: jbrondex

Description:
------------
As Script_PlotProfile_DamageChi_VariousParameters.py, this file makes subplots at different times and for different values of sigma_threshold of damage profiles obtained at
a given forage of the Taconnaz glacier. Damage profiles are represented for different set of parameters (Sigmath, B, lambdah).
The background of each subplot as a specific color depending on wether chi is positive or negative. This way the band of
positive chi (positive source of damage) is clearly highlighted.
In addition, there are 3 sub-subplots for each subplot: one for sigmath constant, one for sigmath as a linear function of Relative Density
(Sigmath_ice*RelDens), and one for the non-linear parameterization of sigmath as a function of Relative Density.
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
from scipy import interpolate

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))

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
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

fig = plt.figure(figsize=(40, 20))
fig.text(0.5, 0.04, r'Damage', ha='center', fontsize=34)
fig.text(0.075, 0.5, r'Depth (m)', va='center', rotation='vertical', fontsize=34)
outer = gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2) ##9 Main Subplots, each divided in 2

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    ##pathroot_mycode = Path('/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLag_DamageSemiLag/ForageOutput/.')
    pathroot_mycode = Path('/home/brondexj/Documents/Taconnaz/MyTaconnaz/MyTaco_Porous_DensSemiLag_DamageSemiLag/ForageOutput/.')
    ### Number of Partitions
    Npartition = 4
    ### Number of the forage of interest (FROM 1 TO 11):
    ForageNumber = 11
    ### Limits of Damage at forage of interest for xlabel limits
    Damage_min = 0
    Damage_max = 0.05
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'Year', 'ForageNumber', 'NodeNumber', 'X', 'Y', 'Z', 'Depth', 'Damage', 'Chi', 'DensRel']
    ###List all combinations of Method and TimeStep to plot on a subplot (with corresponding line style/color)
    B_List = ['20', '20', '10', '10']
    Lambdah_List = ['00','05','00','05']
    Color_List = [RedStar, VeryDarkBrown, RedStar, VeryDarkBrown]
    LineStyle_List = [':', ':', '-', '-']
    LineSize_List = [4, 4, 3, 3]
    Sigmath_List = ['005', '001', '0005']
    SigmathLaw_List = ['Cst', 'FuncRelDens', 'FuncRelDens_NonLine']
    ### Years of interest (one subplot per year):
    Year_List =[1, 10 ,19]
    ### GO FOR IT:
    ### Give the name of Forage Number to the Figure
    Name_Fig = 'Forage N°{}'.format(ForageNumber)
    fig.suptitle(Name_Fig, fontsize=36, weight='bold')
    ### Do for first Sigmath = cst and then Sigmath = f(RelDens) -> Each one of the two subsubplot per subplot
    for j,SigmathLaw in enumerate(SigmathLaw_List):
        #### Do for each value of sigmath (rows of main subplots)
        for i,Sigmath in enumerate(Sigmath_List):
            #### Plot current combination for each of the considered years (columns of main subplots)
            for k,year in enumerate(Year_List):
                l = 3*i + k ###Index of outer subplot
                ### Get proper subplot
                axes = np.empty(shape=(1, 3), dtype=object)  ##This will be for subsubplots
                inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[l], hspace=0, wspace=0)
                ax = plt.Subplot(fig, inner[j], sharey=axes[0, 0])
                if k >0:
                    ax.yaxis.set_visible(False)
                if i <2:
                    ax.xaxis.set_visible(False)
                if j > 0:
                    ax.yaxis.set_visible(False)
                if i == 0: ###Add subtitles to subsubplots corresponding to years
                    # Add ghost axes and titles on columns
                    ax_year = fig.add_subplot(inner[:])
                    ax_year.axis('off')
                    ax_year.set_title('@{}a'.format(year), fontsize=26, weight='bold')
                if k == 0:
                    if i == 0:
                        ax.set_ylabel(r'$\sigma_{th,ice} = 0.05 \mathrm{MPa}$', fontsize=24)
                    elif i==1:
                        ax.set_ylabel(r'$\sigma_{th,ice} = 0.01 \mathrm{MPa}$', fontsize=24)
                    else:
                        ax.set_ylabel(r'$\sigma_{th,ice} = 0.005 \mathrm{MPa}$', fontsize=24)
                    ax.yaxis.set_label_coords(-0.3, 0.5)


                # if j < 2:
                #     ax.set_xlim([Damage_min, Damage_max])
                # elif j == 2:
                #     ax.set_xlim([0.3, 1])
                ax.set_ylim([Damage_min, Damage_max])
                ax.tick_params(labelsize=24)  # fontsize of the tick labels
                ##   ax.grid(True)
                ax.yaxis.set_major_formatter(formatter)
                axes[0, j] = ax
 ##             if i == 0:
 ##                 ax.set_title('@{}a'.format(year), fontsize=34, weight='bold')
                year = year + 1
                ### Load data for all combination in a DataFrame
                for (B, Lambdah, LineCol, LineSize, LineSty) in zip(B_List, Lambdah_List, Color_List, LineSize_List, LineStyle_List):
                    ### Open file corresponding to each partition and then concatenate in single DataFrame
                    Data = []
                    for part in range(Npartition):
                        file_name = 'forage_MyTaco_Porous_SLDensDam_FromStdyD_tsp01_Out1a_20a_Sigth{}_B{}_Sigth{}_Lambdah{}_.dat.{}'.format(
                            SigmathLaw, B, Sigmath, Lambdah, part)
                        file = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names, delim_whitespace=True)
                        Data.append(file)
                    Data = pd.concat(Data)
                    #### Keep only data for the forage of interest and in the ice body (DensRel > 0)
                    Data_Forage = Data[(Data['ForageNumber'] == ForageNumber) & (Data['DensRel'] > 0)]
                    if Data_Forage[Data_Forage['Year']==year]['Damage'].empty:  ####If no results for the considered year (simulation still running) then go to next case
                        continue
                    print('Version:', SigmathLaw, 'Sigmath = ', Sigmath, ',B =', B, ', Lambdah =', Lambdah, ', Time Step n° ',Data_Forage[Data_Forage['Year']==year]['Timestep'].values[0])
                    # Get Variables
                    Depth = Data_Forage[Data_Forage['Year'] == year]['Depth'].values
                    Damage = Data_Forage[Data_Forage['Year'] == year]['Damage'].values
                    Chi = Data_Forage[Data_Forage['Year'] == year]['Chi'].values
                    DensRel = Data_Forage[Data_Forage['Year'] == year]['DensRel'].values
                    #### Plot background using chi values
                    #Reinterpolate Chi on refined Depth mesh
                    Chi_interp_func=interpolate.interp1d(-Depth, Chi,kind='linear') ##linear function to inperpolate chi in between mesh nodes
                    Depth_Refined=np.linspace(0, Depth.max(),1000) ###Refined 1D grid
                    Chi_Refined=Chi_interp_func(-Depth_Refined) ###Chi Reinterpolated on finer 1D mesh
                    # define a binary colormap (one color for negative chi, one color for positive chi)
                    cmap = colors.ListedColormap([BlueGreyBackGround, LightBrownBackGround])
                    # create a normalize object the describes the limits of
                    # each color
                    bounds = [-1, 0.0, 1]
                    norm = colors.BoundaryNorm(bounds, cmap.N)
                    #### Now plot damage profiles:
                    ax.set_ylim([np.min(- Depth) , 0])
                    # if j < 2:
                    #     ax.plot(Damage,- Depth, color=LineCol, linestyle=LineSty, linewidth=LineSize)
                    #     c = ax.pcolorfast([Damage_min, Damage_max],-Depth_Refined, Chi_Refined[np.newaxis].T, cmap=cmap, norm=norm)
                    # elif j == 2:
                    #     ax.plot(DensRel, - Depth, color='r', linestyle='-', linewidth=LineSize)
                    ax.plot(Damage, - Depth, color=LineCol, linestyle=LineSty, linewidth=LineSize)
                    c = ax.pcolorfast([Damage_min, Damage_max],-Depth_Refined, Chi_Refined[np.newaxis].T, cmap=cmap, norm=norm)
                    fig.add_subplot(ax)
    #### DUMMY plot for legend
    B_List_Name = ['2', '2', '1', '1']
    Lambdah_List_Name = ['0.0','0.5','0.0','0.5']
    for i in range(np.size(B_List)):
        ax.plot(np.NaN, np.NaN, label=r'B = {}, $\lambda_h$ = {}'.format(B_List_Name[i], Lambdah_List_Name[i]), color=Color_List[i], linewidth=4, linestyle=LineStyle_List[i])
    fig.subplots_adjust(bottom=0.15, wspace=0.25)
    fig.legend(loc='lower center', fancybox=True, shadow=True,  fontsize=30, ncol=2)


    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'DamageProfile_Forage{}_WithFuncRelDens_NonLine'.format(ForageNumber)
    # name_path = '/home/brondexj/Documents/Taconnaz/MyTaconnaz/MyTaco_Porous_DensSemiLag_DamageSemiLag/Figures/.'
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()