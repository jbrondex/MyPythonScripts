"""
@author: jbrondex

Description:
------------
This scripts is used to produce the result files summarizing main results of the teady simulation regarding the tunnel/channel to drain the Grande Motte glacier lake.
"""

################################################################################
# import #####
################################################################################
from pathlib import Path
import numpy as np
import pandas as pd  ###To treat the csv files
import csv ###To write results in a csv file

####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN PART OF THE CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
if __name__ == "__main__":
    #### USER DEFINED PARAMETERS :####
    ### Where to find the output files:
    pathroot_mycode = Path('/home/brondexj/BETTIK/GrandeMotteTignes/.')
    #########################################################################
    ####   BELOW CODE LINES TO PRODUCE RESULTS FILE FOR STEADY SIMUS     ####
    #########################################################################
    ##BE AWARE: the eigen stress are not in the same order for elastic and viscous (eigenstress 1 and 3 for the viscous are resp. min and max)
    Col_Names_Steady = ['Tsp','BC','Node','X','Y','Z','Vx','Vy','SigmaI','Vn','Sigma_nn']
    ###List all combinations of transects, shapes and widths
    Transect_List = ['Transect1', 'Transect2', 'Transect3', 'Transect4', 'Transect5', 'Transect6']
    Shape_List =['Rectangle', 'Circle', 'Ovoide']
    Width_List = ['W100cm','W150cm', 'W200cm']

    ###Create a file that summarizes the results of the steady simulations
    filename_results_steady = './RunSteady_Tunnel/Results_SimuSteady.csv'
    ###Create columns in the results file ('Transect', 'Shape', 'Width', 'SigmaI min', 'SigmaI max', SigmaI mean', 'Vn min', 'Vn max', Vn mean')
    Header_Steady = ['Transect', 'Shape', 'Width', 'SigmaI_max [kPa]', 'SigmaI_mean [kPa]', 'Vn_max [cm/d]']
    with open(pathroot_mycode.joinpath(filename_results_steady), 'w', newline ='') as ResultSteady:
        write = csv.writer(ResultSteady, delimiter=" ")
        write.writerow(Header_Steady)
        ### Start the loop on transects
        for Transect in Transect_List:
            ### Start the loop on shapes
            for Shape in Shape_List:
                if Shape == 'Ovoide':
                    ####Name of the elmer output file
                    filename_output = 'RunSteady_Tunnel/RunSteady_{0}_{1}/ScalarOutput/TunnelOutput_RunSteady_Tunnel_{0}_{1}_TunnelFreeSurf_.dat'.format(Transect, Shape)
                    df = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Steady, delim_whitespace=True)
                    SigmaI_min = np.min(df['SigmaI'])*1000
                    SigmaI_max = np.max(df['SigmaI'])*1000
                    SigmaI_mean = np.mean(df['SigmaI'])*1000
                    Vn_min = np.min(df['Vn'])*100/365.25
                    Vn_max = np.max(df['Vn'])*100/365.25
                    Vn_mean = np.mean(df['Vn'])*100/365.25
                    print(Transect, ' ', Shape, ' ', 'SigmaImax =', SigmaI_max, ' SigmaImean =', SigmaI_mean )
                    Result_Line=[Transect, Shape, '-', SigmaI_max, SigmaI_mean, Vn_max]
                    write.writerow(Result_Line)
                else:
                    for Width in Width_List:
                        ####Name of the elmer output file
                        filename_output = 'RunSteady_Tunnel/RunSteady_{0}_{1}_{2}/ScalarOutput/TunnelOutput_RunSteady_Tunnel_{0}_{1}_{2}_TunnelFreeSurf_.dat'.format(Transect, Shape, Width)
                        df = pd.read_csv(pathroot_mycode.joinpath(filename_output), names=Col_Names_Steady, delim_whitespace=True)
                        SigmaI_min = np.min(df['SigmaI']) * 1000
                        SigmaI_max = np.max(df['SigmaI']) * 1000
                        SigmaI_mean = np.mean(df['SigmaI']) * 1000
                        Vn_min = np.min(df['Vn']) * 100 / 365.25
                        Vn_max = np.max(df['Vn']) * 100 / 365.25
                        Vn_mean = np.mean(df['Vn']) * 100 / 365.25
                        print(Transect, ' ', Shape, ' ', Width, ' SigmaImax =', SigmaI_max, ' SigmaImean =', SigmaI_mean )
                        Result_Line = [Transect, Shape, Width, SigmaI_max, SigmaI_mean, Vn_max]
                        write.writerow(Result_Line)



