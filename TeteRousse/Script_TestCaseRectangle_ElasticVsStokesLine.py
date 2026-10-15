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

fig3, axes2 = plt.subplots(2, 3, figsize=(60, 60))
fig3.subplots_adjust(hspace=0.2,wspace=0.35)
for i in range(0,2):
    for j in range(0,3):
        ax = axes2[i,j]
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

fig4, axes4 = plt.subplots(2, 3, figsize=(60, 60))
fig4.subplots_adjust(hspace=0.2,wspace=0.35)
for i in range(0,2):
    for j in range(0,3):
        ax = axes4[i,j]
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
    pathroot_mycode = Path('/home/brondexj/Documents/TestsBenchmark_Elmer/TestCase_LinearElasticVSStokesLine/ScalarOutput/.')
    ##BE AWARE: the eigen stress are not in the same order for elastic and viscous (eigenstress 1 and 3 for the viscous are resp. min and max)
    Col_Names_Viscous = ['Tsp','Line','Node','X','Y','Z','SigmaMin','SigmaMiddle','SigmaMax']
    ### For the elastic case
    ###List all combinations of parameters in the case of the elastic solver
    ##BE AWARE: the eigen stress are not in the same order for elastic and viscous (Principal Stress 1 and 3 for the elastic are resp. max and min)
    Col_Names_Elastic = ['Tsp','Line','Node','X','Y','Z','SigmaMax','SigmaMiddle','SigmaMin']
    PoissonName_List = ['030', '035', '040', '045', '049', '04999']
    Poisson_List =[0.3,0.35,0.4,0.45,0.49,0.4999]
    Case_List = ['','_CompressionTest']
    Young_List = ['1', '9']
    LineStyle_List = [':', '-']

    ### Start the loop on the horizontal (Line 1) and Vertical (Line 2) profiles (= rows of subplot)
    for Case in Case_List:
        ### Load the viscous case (reference !!)
        file_name = 'Output_Test_StokesLine_BCBottomLeftNode{}_.dat'.format(Case)
        df_Viscous = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names_Viscous, delim_whitespace=True)
        ###Same for non-linear case
        file_name_NL = 'Output_Test_StokesNonLine_BCBottomLeftNode{}_.dat'.format(Case)
        df_Viscous_NL = pd.read_csv(pathroot_mycode.joinpath(file_name_NL), names=Col_Names_Viscous,delim_whitespace=True)
        for i,Line in enumerate([1,2]):
            ### We select only outputs along the current line
            df_Viscous_tmp=df_Viscous[df_Viscous['Line']==Line]
            df_Viscous_NL_tmp=df_Viscous_NL[df_Viscous_NL['Line']==Line]
            ## Loop on the 3 principle stresses (columns of subplot)
            for j,Stress in enumerate(['SigmaMax', 'SigmaMiddle', 'SigmaMin']):
                ##Get proper ax
                if Case == '':
                    ax=axes[i,j]
                elif Case == '_CompressionTest':
                    ax=axes2[i,j]
                    if i == 0:
                        ax.set_ylim([-1.1, 0.1])
                    elif i == 1:
                        ax.set_xlim([-1.1, 0.1])
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
                        filename = 'Output_Test_Elastic_BCBottomLeftNode{}_Poisson{}_Young{}GPa_.dat'.format(Case,PoissonName, Young)
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
    fig3.colorbar(cmx.ScalarMappable(norm=cNorm, cmap=cm), ax=axes2[0,:], ticks=[0.3,0.35,0.4,0.45,0.5],orientation='vertical', label=r'Poisson Coeff')
    fig3.colorbar(cmx.ScalarMappable(norm=cNorm, cmap=cm), ax=axes2[1,:], ticks=[0.3,0.35,0.4,0.45,0.5],orientation='vertical', label=r'Poisson Coeff')

    #### FROM THERE FIGURE 4 : STRESS FIELD FOR STOKES LINE WITH FRICTION AND VARIOUS A  ####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/Documents/TestsBenchmark_Elmer/TestCase_StokesLineVariousA_WithFriction/ScalarOutput/.')
    ###List of arrhenius factor tested
    A_List = ['02', '04', '08', '12']
    LineStyle_List = [':', '-.', '--', '-']
    Col_List = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    SlipCoeff = 'Cw01' ##'NoSlid'##'Cw02'
    ### Start the loop on the horizontal (Line 1) and Vertical (Line 2) profiles (= rows of subplot)
    for (A, Sty, Col) in zip(A_List,LineStyle_List,Col_List):
        ### Load the viscous case (reference !!)
        file_name = 'Output_Test_StokesLine_{}_A{}.dat'.format(SlipCoeff,A)
        df_Viscous = pd.read_csv(pathroot_mycode.joinpath(file_name), names=Col_Names_Viscous, delim_whitespace=True)
        for i, Line in enumerate([1, 2]):
            ### We select only outputs along the current line
            df_Viscous_tmp = df_Viscous[df_Viscous['Line'] == Line]
            ## Loop on the 3 principle stresses (columns of subplot)
            for j, Stress in enumerate(['SigmaMax', 'SigmaMiddle', 'SigmaMin']):
                ##Get proper ax
                ax = axes4[i, j]
                ##Add title to columns
                if i == 0:
                    if j == 0:
                        ax.set_title(r'$\sigma_I$', fontsize=21, weight='bold')
                    elif j == 1:
                        ax.set_title(r'$\sigma_{II}$', fontsize=21, weight='bold')
                    else:
                        ax.set_title(r'$\sigma_{III}$', fontsize=21, weight='bold')
                ##Plot the reference stress profile (i.e. viscous stress) on top
                if i == 0:
                    df_Viscous_tmp.sort_values(by=['X'], inplace=True)
                    ax.plot(df_Viscous_tmp['X'].values, df_Viscous_tmp[Stress].values, color=Col, linewidth=2.2,linestyle=Sty)
                else:
                    df_Viscous_tmp.sort_values(by=['Y'], inplace=True)
                    ax.plot(df_Viscous_tmp[Stress].values, df_Viscous_tmp['Y'].values, color=Col, linewidth=2.2,linestyle=Sty)


    plt.show()






