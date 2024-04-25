"""
@author: jbrondex

Description:
------------
This file makes 3 plots to explore the parameterization of sigma_threshold as a function of relative density in firn.
First plot shows the evolution of Sigmath_firn as a function of relative density for three value of sigmath_ice. We
explore three parameterization: constant with Sigmath_firn = sigmath_ice, linear with Sigmath_firn = RelDens * sigmath_ice,
non-linear with Sigmath_firn = Factor * sigmath_ice and Factor=f(a,b,n)
Second plot is actually made of 11 subplots (one per forage of Taconnaz glacier). It shows profiles of Sigmath_firn for the three values
of sigmath_ice and for the three considered parameterization as well as the distribution of SigmaI (maximum pronciple stress) as obtained
after 20 years of simu in Taconnaz glacier.
Third plot is a zoom-in of second plot in the near-surface region.
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

###Fig 1 is Sigmath as function of RelDens for different parameterization
fig1 = plt.figure(1, figsize=(40, 20))
plt.ylabel(r'$\sigma_{th}$ (MPa)', fontsize=34)
plt.xlabel(r'Rel Dens', fontsize=34)
plt.xlim([0.0, 1.1])
plt.tick_params(labelsize=24)  # fontsize of the tick labels
plt.grid(True)


###Fig 2 is profile of Sigmath as function of depth for each of the 11 forage
rows=3
col=4
# ###Parameters for save figure
# fig2, axes = plt.subplots(rows, col, figsize=(80, 60), sharex=True, sharey=True)
# fig2.text(0.5, 0.07, r'$\sigma_{th}$ (MPa)', fontsize=90, ha='center')
# fig2.text(0.08, 0.5, r'Depth (m)', fontsize=90, va='center', rotation='vertical')
# fig2.subplots_adjust(hspace=0.1,wspace=0.12)
# axes[-1, -1].axis('off') ##Because only 11 forage and not 12
# for i in range(0,rows):
#     for j in range(0,col):
#         ax = axes[i,j]
#         ax.set_xlim([0.0, 0.2])
#         ax.tick_params(labelsize=58)  # fontsize of the tick labels
#         ax.grid(True)
#         ax.yaxis.set_major_formatter(formatter)

###Parameters if figure not saved
fig2, axes = plt.subplots(rows, col, figsize=(40, 30), sharex=True, sharey=True)
fig2.text(0.5, 0.05, r'$\sigma_{th}$ (MPa)', fontsize=16, ha='center')
fig2.text(0.08, 0.5, r'Depth (m)', fontsize=16, va='center', rotation='vertical')
fig2.subplots_adjust(hspace=0.2,wspace=0.2)
axes[-1, -1].axis('off') ##Because only 11 forage and not 12
for i in range(0,rows):
    for j in range(0,col):
        ax = axes[i,j]
        ax.set_xlim([0.0, 0.2])
        ax.tick_params(labelsize=12)  # fontsize of the tick labels
        ax.grid(True)
        ax.yaxis.set_major_formatter(formatter)

###Parameters if figure not saved (Zoom-In)
fig3, axes3 = plt.subplots(rows, col, figsize=(40, 30), sharex=True, sharey=True)
fig3.text(0.5, 0.05, r'$\sigma_{th}$ (MPa)', fontsize=16, ha='center')
fig3.text(0.08, 0.5, r'Depth (m)', fontsize=16, va='center', rotation='vertical')
fig3.subplots_adjust(hspace=0.2,wspace=0.2)
axes3[-1, -1].axis('off') ##Because only 11 forage and not 12
for i in range(0,rows):
    for j in range(0,col):
        ax = axes3[i,j]
        ax.set_xlim([0.0, 0.005])
        ax.set_ylim([-5.0, 0.0])
        ax.tick_params(labelsize=12)  # fontsize of the tick labels
        ax.grid(True)
        ax.yaxis.set_major_formatter(formatter)
###############################
# Define function a and b #####
###############################
def FuncA(RelDens):
    if RelDens > 0.99:
        a=1.0
    elif (RelDens <= 0.81):
        dd= RelDens
        if (RelDens < 0.4):
            dd = 0.4
        a=np.exp(13.22240 - 15.78652 * dd)
    else:
        a= 1.0 + (2/3) * (1 -RelDens)
        a = a / (RelDens**1.5)
    return a

def FuncB(RelDens):
    if RelDens > 0.99:
        b=0.0
    elif (RelDens <= 0.81):
        dd= RelDens
        if (RelDens < 0.4):
            dd = 0.4
        b=np.exp(15.09371 - 20.46489 * dd)
    else:
        b= (1.0 - RelDens)**(1/3)
        b=(3/4) * (b/ (3 * (1 -b)))**1.5
    return b

def Sigrist2006_SigmaTh(RelDens):
    Th = 0.24* RelDens**2.44
    return Th

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ### Fig1: Simply SigmathPorous as a function of RelDens for param line and param rheo for different sigmath_ice
    n = 3 ###Glen Parameter
    SigmathIce_List = [0.2, 0.1, 0.05]
    LineStyle_List = ['-', '--', ':']
    # LineSize_List = [7, 7, 7] ##Parameters if figure saved
    LineSize_List = [3, 3, 3]
    RelDens_Dummy = np.linspace(0,1.1,100)
    for (SigmathIce, LineStyle, LineSize) in zip(SigmathIce_List, LineStyle_List, LineSize_List):
        Sigmath_Line=RelDens_Dummy*SigmathIce
        Sigmath_Rheo = []
        for k in range(len(RelDens_Dummy)):
            value=(FuncA(RelDens_Dummy[k])+FuncB(RelDens_Dummy[k])/3)**(-((n+1)/(2*n)))*SigmathIce
            Sigmath_Rheo.append(value)
        ### For Parameterization of Sigrist 2006 (only valid for snow dens between 80 and 350
        RelDens_Dummy_Filtered=RelDens_Dummy[RelDens_Dummy[:] < 350/917]
        RelDens_Dummy_Filtered=RelDens_Dummy_Filtered[RelDens_Dummy_Filtered[:]>80/917]
        Sigmath_Sigrist = []
        for k in range(len(RelDens_Dummy_Filtered)):
            value=Sigrist2006_SigmaTh(RelDens_Dummy_Filtered[k])
            Sigmath_Sigrist.append(value)
        fig1 = plt.figure(1)
        plt.axhline(y=SigmathIce, color='k', linestyle=LineStyle, linewidth=2)
        plt.plot(RelDens_Dummy,Sigmath_Line,color='r',linestyle=LineStyle, linewidth=LineSize)
        plt.plot(RelDens_Dummy, Sigmath_Rheo, color='darkblue',linestyle=LineStyle, linewidth=LineSize)
        plt.plot(RelDens_Dummy_Filtered, Sigmath_Sigrist, color=LightBrown,linestyle=LineStyle, linewidth=LineSize)
    # plt.show()

    ###Fig2: for each forage, sigmath as function of depth
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   SUBPLOTS f(Depth) FOR ALL FORAGE LOAD ALL DATA AND PLOT    ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    ### Where to find the output files:
    ##pathroot_mycode = Path('/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLag_DamageSemiLag/ForageOutput/.')
    pathroot_mycode = Path('/home/brondexj/Documents/Taconnaz/MyTaconnaz/MyTaco_Porous_DensSemiLag_DamageSemiLag/ForageOutput/.')
    ### Number of Partitions
    Npartition = 4
    ### Number of the forage of interest (FROM 1 TO 11):
    ForageNumber = 11
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'Year', 'ForageNumber', 'NodeNumber', 'X', 'Y', 'Z', 'Depth', 'Damage', 'Chi', 'DensRel', 'SigmaI']
    Name_Fig = 'Forage N°{}'.format(ForageNumber)

    ### Open file with dens profiles corresponding to each partition and then concatenate in single DataFrame
    Data = []
    for part in range(Npartition):
        file_name = 'forage_MyTaco_Porous_SLDensDam_FromStdyD_tsp01_Out1a_20a_SigthCst_B10_Sigth01_Lambdah00_PrincipStressOut_.dat.{}'.format(part)
        file = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names, delim_whitespace=True)
        Data.append(file)
    Data = pd.concat(Data)
    ###Keep only last year of simu (1 because I do a restart for a single tsp to return fields that I forgot)
    Data = Data[(Data['Year'] == 1)]
    ### Do one subplot per forage
    for i in range(rows): ###Lines of subplot
        for j in range(col): ###col of subplot
            if i==rows-1 and j==col-1:  ###skip last subplot as there are 11 forage, not 12
                continue
            Forage= 4*i+1+j ##Les forages sont numérotés en partant de 1
            ### Get proper axe of subplot
            ax = axes[i,j]
            ax.set_title('Forage N°{}'.format(Forage), fontsize =14) ## if save figure then fontsize=60
            ### Same for zoom-in
            ax3 = axes3[i,j]
            ax3.set_title('Forage N°{}'.format(Forage), fontsize =14) ## if save figure then fontsize=60
            #### Keep only data for the forage of interest and in the ice body (DensRel > 0)
            Data_Forage = Data[(Data['ForageNumber'] == Forage) & (Data['DensRel'] > 0)]
            # Get Variables and calculates corresponding Sigmath
            DensRel = Data_Forage['DensRel'].values
            Depth= Data_Forage['Depth'].values
            SigmaI = Data_Forage['SigmaI'].values
            for (SigmathIce, LineStyle, LineSize) in zip(SigmathIce_List, LineStyle_List, LineSize_List):
                Sigmath_Line = DensRel * SigmathIce
                Sigmath_Rheo = []
                for k in range(len(DensRel)):
                    value = (FuncA(DensRel[k]) + FuncB(DensRel[k]) / 3) ** (-((n + 1) / (2 * n))) * SigmathIce
                    Sigmath_Rheo.append(value)
                ax.axvline(x=SigmathIce, color='k', linestyle=LineStyle, linewidth =1.5) ##If figure saved linewidth=5
                ax.plot(Sigmath_Line, -Depth, color='r', linestyle=LineStyle, linewidth=LineSize)
                ax.plot(Sigmath_Rheo, -Depth, color='darkblue', linestyle=LineStyle, linewidth=LineSize)
                ax.plot(SigmaI, -Depth, color='forestgreen', linestyle=LineStyle, linewidth=LineSize)
                ###Same for zoom-in
                ax3.axvline(x=SigmathIce, color='k', linestyle=LineStyle, linewidth =1.5) ##If figure saved linewidth=5
                ax3.plot(Sigmath_Line, -Depth, color='r', linestyle=LineStyle, linewidth=LineSize)
                ax3.plot(Sigmath_Rheo, -Depth, color='darkblue', linestyle=LineStyle, linewidth=LineSize)
                ax3.plot(SigmaI, -Depth, color='forestgreen', linestyle=LineStyle, linewidth=LineSize)


    #### DUMMY plot for legend
    ax.plot(np.NaN, np.NaN, label=r'$\sigma_I$', color='forestgreen', linewidth=3, linestyle='-')
    ax.plot(np.NaN, np.NaN, label=r'$\sigma_{th}$ Linear', color='r', linewidth=3, linestyle='-')
    ax.plot(np.NaN, np.NaN, label=r'$\sigma_{th}$ Non Linear', color='darkblue', linewidth=3, linestyle='-')
    for i in range(np.size(SigmathIce_List)):
        thice='{th,ice}'
        ax.plot(np.NaN, np.NaN, label=r'$\sigma_{}$  = {} MPa'.format(thice,str(SigmathIce_List[i])), color='k', linewidth=3, linestyle=LineStyle_List[i])
   ## fig2.legend(loc='lower left', bbox_to_anchor=(0.71,0.18), fancybox=True, shadow=True,  fontsize=58, ncol=2)
    fig2.legend(loc='lower left', bbox_to_anchor=(0.73,0.10), fancybox=True, shadow=True,  fontsize=16, ncol=1)
    ax3.plot(np.NaN, np.NaN, label=r'$\sigma_I$', color='forestgreen', linewidth=3, linestyle='-')
    ax3.plot(np.NaN, np.NaN, label=r'$\sigma_{th}$ Linear', color='r', linewidth=3, linestyle='-')
    ax3.plot(np.NaN, np.NaN, label=r'$\sigma_{th}$ Non Linear', color='darkblue', linewidth=3, linestyle='-')
    for i in range(np.size(SigmathIce_List)):
        thice='{th,ice}'
        ax3.plot(np.NaN, np.NaN, label=r'$\sigma_{}$  = {} MPa'.format(thice,str(SigmathIce_List[i])), color='k', linewidth=3, linestyle=LineStyle_List[i])
    fig3.legend(loc='lower left', bbox_to_anchor=(0.73,0.10), fancybox=True, shadow=True,  fontsize=16, ncol=1)
    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'SigmathProfile_ForAllForage_'
    # name_path = '/home/brondexj/Documents/Taconnaz/MyTaconnaz/MyTaco_Porous_DensSemiLag_DamageSemiLag/Figures/.'
    # path_output_fig = Path(name_path)
    # fig2.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()