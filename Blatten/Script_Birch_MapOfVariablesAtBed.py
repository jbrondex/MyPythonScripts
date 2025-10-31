"""
@author: jbrondex

Description:
------------
This file aims at plotting figures corresponding to map of relevant variable at bed at given time
and a vertical profile along y = y_borehole 2

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch
from matplotlib.colors import LogNorm

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib.ticker import FormatStrFormatter
from matplotlib.dates import date2num
import matplotlib.gridspec as gridspec
import rasterio
from rasterio.windows import from_bounds

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
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y)
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

def Interpolate_field_slice(df, field_name, x, z): ###Returned field interpolated from mesh grid to point (x,z) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    zi = df['Z'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, zi), field, (x, z), rescale=True)
    return result

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
    VariableToPlot = 'Sigma_nn' ###What variable to plot (to make script faster): 'Velocity', 'Sigma_nn', 'Frottement', 'Cavitation'
    ###Parameter of Simu
    C = 0.3
    ### velocity max for plot of ub (to adapt depending on C)
    velocitymax = 30  ##in m/d !
    #   Timestep_size = 0.01/365 ### In years !!
    Timestep_size = 0.01 / 365  ### In years !!
    OutputInterval = 10  ###Execution interval of the solver producing the simu output
    OutputOfPlot = 3 ##[0, 1, 3, 7]  ###At which output (count) do we want to make the plot (correspond to the day of simu)
    ###Coordinate of bed point at which we want to plot the variables vs time (this is done in the other python script)
    xplot = [2630580, 2630620.9, 2630653.6]
    yplot = [1139265, 1139221.7, 1139188.3]
    colplot = ['#00FF66', '#8D5FD3', '#FFD42A']
    Stars = False  ##False ## Do we want to plot stars ?
    ### Size of the stars on the plot
    Marker_Size = 18  ### Size of markers to show where vertical profiles of variables vs time are plotted
    ####Plot parameters###
    ### Contours we want to plot for debris thickness
    DebrisContour_Levels = [10, 20, 30, 40, 50, 60]
    DebrisContour_LabelPositions = [(2630.73, 1139.097), (2630.587, 1139.089), (2630.877, 1139.018),
                                    (2630.917, 1138.978), (2630.736, 1138.978), (2630.853, 1138.930)]
    DebrisContour_Color = 'lightcoral'
    DebrisContour_Thick = 0.9
    ###Defines bounds of the plot
    xmin, xmax = 2630250, 2631250
    ymin, ymax = 1138750, 1139500
    ###Open Data corresponding to grounded mask
    # Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
    # str_Cvalue = f"{int(C * 10):02d}"
    # Filename = 'MyBirch_RS_C{}_WithD_Cont__bed.dat'.format(str_Cvalue)
    Pathroot = Path(
        '/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/4_RateAndState_WithDamage_WithContact/GlacierOut/.')
    str_Cvalue = f"{int(C * 10):02d}"
    Filename = 'MyBirch_RS_C{}_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard__bed.dat'.format(str_Cvalue)
    Col_Names = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'H', 'N', 'theta_dagger', 'theta', 'taub',
                 'nodearea', 'fx', 'fy', 'fz', 'u', 'v', 'w', 'nx', 'ny', 'nz', 'Chi', 'GM']

    ### START OF SCRIPT###
    ###Load glacier contour
    Pathroot_Contour = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/Visu_Bed/.')
    Filename_GlacierContour = 'ContourGlacier_FromMesh3D.dat'
    Df_GlacierContour = pd.read_csv(Pathroot_Contour.joinpath(Filename_GlacierContour), names=['X', 'Y'], sep=r'\s+')
    xc = Df_GlacierContour['X'].values
    xc=np.append(xc,xc[0]) ##step required to close the contour
    yc = Df_GlacierContour['Y'].values
    yc = np.append(yc, yc[0])  ##step required to close the contour

    ### Load Orthophoto for Background
    # --- Load orthophoto
    Pathroot_Ortho = Path('/home/brondexj/Documents/Expertises/Blatten/OrthoPhotos/.')
    Filename_Ortho = 'OrthoPostCollapse.tif'
    with rasterio.open(Pathroot_Ortho.joinpath(Filename_Ortho)) as src:
        window = from_bounds(xmin, ymin, xmax, ymax, src.transform)
        ortho_crop = src.read(1, window=window)  # for grayscale
        # ortho_crop = src.read([1, 2, 3], window=window) # for RGB
        transform = src.window_transform(window)
    # --- Compute extent in km
    extent_km = [xmin / 1000, xmax / 1000, ymin / 1000, ymax / 1000]

    ### Load DEM of surface of glacier and rock deposit to build contour line of rock debris thickness
    Pathroot_DEMs= Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/1_Mesh/MESH2025/.')
    Filename_SurfGlacierDEM = 'DEMSurfPreIce.dat'
    Df_SurfGlacierDEM = pd.read_csv(Pathroot_DEMs.joinpath(Filename_SurfGlacierDEM), names=['X', 'Y', 'Z'], sep=r'\s+')
    Filename_SurfDebrisDEM = 'DEMSurfPre.dat'
    Df_SurfDebrisDEM = pd.read_csv(Pathroot_DEMs.joinpath(Filename_SurfDebrisDEM), names=['X', 'Y', 'Z'], sep=r'\s+')

    #########################################
    ####  OPEN OUTPUT OF THE SIMULATIONS:####
    #########################################
    ###load file as dataframe
    Df_wholebed_full = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Names, sep=r'\s+')
    #####For the with contact case, several restart have been required: build a single dataframe from all restart
    ###number of restart required depends on value of C
    ###number of restart required depends on value of C
    if C == 0.7:
        Filename_Rstt1 = 'MyBirch_RS_C07_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp2901__bed.dat'
        Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt1], ignore_index=True)
    elif C == 0.6:
        ###For this case, there is overlap between end of simu R0 and beginning of simu R1
        ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 2501)
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 2501].reset_index(drop=True)
        Filename_Rstt1 = ('MyBirch_RS_C06_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp2501__bed.dat')
        Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt1], ignore_index=True)
        ###Same thing for restart R2 (tsp 1000 of R1 so 3501 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 3501].reset_index(drop=True)
        Filename_Rstt2 = 'MyBirch_RS_C06_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp1000__bed.dat'
        Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt2], ignore_index=True)

    elif C == 0.5:
        ###For this case, there is overlap between end of simu R0 and beginning of simu R1
        ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 1601)
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 1601].reset_index(drop=True)
        Filename_Rstt1 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp1601__bed.dat'
        Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt1], ignore_index=True)
        ###Same thing for restart R2 (tsp 600 of R1 so 2201 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 2201].reset_index(drop=True)
        Filename_Rstt2 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp600__bed.dat'
        Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt2], ignore_index=True)
        ###Same thing for restart R3 (tsp 500 of R2 so 2701 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 2701].reset_index(drop=True)
        Filename_Rstt3 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt3_tsp500__bed.dat'
        Df_wholebed_Rstt3 = pd.read_csv(Pathroot.joinpath(Filename_Rstt3), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt3 = Df_wholebed_Rstt3[Df_wholebed_Rstt3["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt3["Timestep"] = Df_wholebed_Rstt3["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt3["Count"] = Df_wholebed_Rstt3["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt3], ignore_index=True)
        ###Same thing for restart R4 (tsp 700 of R3 so 3401 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 3401].reset_index(drop=True)
        Filename_Rstt4 = 'MyBirch_RS_C05_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt4_tsp700__bed.dat'
        Df_wholebed_Rstt4 = pd.read_csv(Pathroot.joinpath(Filename_Rstt4), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt4 = Df_wholebed_Rstt4[Df_wholebed_Rstt4["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt4["Timestep"] = Df_wholebed_Rstt4["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt4["Count"] = Df_wholebed_Rstt4["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt4], ignore_index=True)
    elif C == 0.4:
        ###For this case, there is overlap between end of simu R0 and beginning of simu R1
        ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 1101)
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 1101].reset_index(drop=True)
        Filename_Rstt1 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp1101__bed.dat'
        Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt1], ignore_index=True)
        ###Same thing for restart R2 (tsp 400 of R1 so 1501 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 1501].reset_index(drop=True)
        Filename_Rstt2 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp400__bed.dat'
        Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt2], ignore_index=True)
        ###Same thing for restart R3 (tsp 600 of R2 so 2101 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 2101].reset_index(drop=True)
        Filename_Rstt3 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt3_tsp600__bed.dat'
        Df_wholebed_Rstt3 = pd.read_csv(Pathroot.joinpath(Filename_Rstt3), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt3 = Df_wholebed_Rstt3[Df_wholebed_Rstt3["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt3["Timestep"] = Df_wholebed_Rstt3["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt3["Count"] = Df_wholebed_Rstt3["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt3], ignore_index=True)
        ###Same thing for restart R4 (tsp 500 of R3 so 2601 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 2601].reset_index(drop=True)
        Filename_Rstt4 = 'MyBirch_RS_C04_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt4_tsp500__bed.dat'
        Df_wholebed_Rstt4 = pd.read_csv(Pathroot.joinpath(Filename_Rstt4), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt4 = Df_wholebed_Rstt4[Df_wholebed_Rstt4["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt4["Timestep"] = Df_wholebed_Rstt4["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt4["Count"] = Df_wholebed_Rstt4["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt4], ignore_index=True)
    elif C == 0.3:
        ###For this case, there is overlap between end of simu R0 and beginning of simu R1
        ###To remove the overlap I remove all lines of simu R0 that are after the restart time of R1 (tsp 701)
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 701].reset_index(drop=True)
        Filename_Rstt1 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt1_tsp701__bed.dat'
        Df_wholebed_Rstt1 = pd.read_csv(Pathroot.joinpath(Filename_Rstt1), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt1 = Df_wholebed_Rstt1[Df_wholebed_Rstt1["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt1["Timestep"] = Df_wholebed_Rstt1["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt1["Count"] = Df_wholebed_Rstt1["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt1], ignore_index=True)
        ###Same thing for restart R2 (tsp 300 of R1 so 1001 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 1001].reset_index(drop=True)
        Filename_Rstt2 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt2_tsp300__bed.dat'
        Df_wholebed_Rstt2 = pd.read_csv(Pathroot.joinpath(Filename_Rstt2), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt2 = Df_wholebed_Rstt2[Df_wholebed_Rstt2["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt2["Timestep"] = Df_wholebed_Rstt2["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt2["Count"] = Df_wholebed_Rstt2["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt2], ignore_index=True)
        ###Same thing for restart R3 (tsp 500 of R2 so 1501 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 1501].reset_index(drop=True)
        Filename_Rstt3 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt3_tsp500__bed.dat'
        Df_wholebed_Rstt3 = pd.read_csv(Pathroot.joinpath(Filename_Rstt3), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt3 = Df_wholebed_Rstt3[Df_wholebed_Rstt3["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt3["Timestep"] = Df_wholebed_Rstt3["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt3["Count"] = Df_wholebed_Rstt3["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt3], ignore_index=True)
        ###Same thing for restart R4 (tsp 400 of R3 so 1901 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 1901].reset_index(drop=True)
        Filename_Rstt4 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt4_tsp400__bed.dat'
        Df_wholebed_Rstt4 = pd.read_csv(Pathroot.joinpath(Filename_Rstt4), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt4 = Df_wholebed_Rstt4[Df_wholebed_Rstt4["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt4["Timestep"] = Df_wholebed_Rstt4["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt4["Count"] = Df_wholebed_Rstt4["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt4], ignore_index=True)
        ###Same thing for restart R5 (tsp 400 of R4 so 2301 from R0) :
        Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Timestep'] <= 2301].reset_index(drop=True)
        Filename_Rstt5 = 'MyBirch_RS_C03_WithD_Cont_SaveScalInt10_tsp001d_ThetaLimited_FullPicard_Rstt5_tsp400__bed.dat'
        Df_wholebed_Rstt5 = pd.read_csv(Pathroot.joinpath(Filename_Rstt5), names=Col_Names, sep=r'\s+')
        # 1. Remove rows where Timestep == 1
        Df_wholebed_Rstt5 = Df_wholebed_Rstt5[Df_wholebed_Rstt5["Timestep"] != 1].copy()
        # 2. Shift timesteps to continue after the last timestep of Df_wholebed
        last_timestep = Df_wholebed_full["Timestep"].max()
        Df_wholebed_Rstt5["Timestep"] = Df_wholebed_Rstt5["Timestep"] + last_timestep - 1
        last_count = Df_wholebed_full["Count"].max()
        Df_wholebed_Rstt5["Count"] = Df_wholebed_Rstt5["Count"] + last_count - 1
        # 3. Concatenate back into Df_wholebed
        Df_wholebed_full = pd.concat([Df_wholebed_full, Df_wholebed_Rstt5], ignore_index=True)

    ###Create a column with time (in years !)
    Df_wholebed_full['Time'] = Timestep_size * Df_wholebed_full['Timestep']
    # We look only the 40 first days : Keep only rows where Time <= 40/365 and reset index
    Df_wholebed_full = Df_wholebed_full[Df_wholebed_full['Time'] <= 40 / 365].reset_index(drop=True)

    ########################################################################
    ###     GET VARIABLES OF INTEREST AT DAY OF INTEREST AND INTERPOLATE ###
    ########################################################################
    print('Dealing with plot for day', OutputOfPlot)
    ### Keep only lines for timestep at which we want to make the plot
    Df_wholebed = Df_wholebed_full[Df_wholebed_full['Count'] == int((OutputOfPlot / (Timestep_size * OutputInterval * 365)) + 1)].copy()
    ### At the moment of this script, there was a bug in the diagnostic of taub in Elmer. Force taub to zero at all floating nodes:
    Df_wholebed.loc[Df_wholebed['GM'] < 0, 'taub'] = 0
    ###Calculate ub from velocity and normal vector
    ###Start by calculating un (should be zero but there is some residual)
    Df_wholebed['un'] = Df_wholebed['u'] * Df_wholebed['nx'] + Df_wholebed['v'] * Df_wholebed['ny'] + Df_wholebed['w'] * Df_wholebed['nz']
    ###Deduce ub from u and un
    Df_wholebed['ub'] = np.sqrt((Df_wholebed['u'] - Df_wholebed['un'] * Df_wholebed['nx']) ** 2 + (Df_wholebed['v'] - Df_wholebed['un'] * Df_wholebed['ny']) ** 2 + (Df_wholebed['w'] - Df_wholebed['un'] * Df_wholebed['nz']) ** 2)
    ### Calculate normal stresses against bedrock by projecting fx,fy,fz on the normal
    Df_wholebed['Sigma_nn'] = (Df_wholebed['fx'] * Df_wholebed['nx'] + Df_wholebed['fy'] * Df_wholebed['ny'] + Df_wholebed['fz'] * Df_wholebed['nz']) / Df_wholebed['nodearea']
    ###Create refined regular grid and interpolate field over grid
    xbed = np.arange(np.floor(np.min(Df_wholebed['X']))-10,np.floor(np.max(Df_wholebed['X']))+10,0.1)
    ybed = np.arange(np.floor(np.min(Df_wholebed['Y']))-10,np.floor(np.max(Df_wholebed['Y']))+10,0.1)
    Xbed, Ybed = np.meshgrid(xbed, ybed)
    ### Interpolate DEMs over same grid (just needed once)
    ZSurfGlacier = Interpolate_field(Df_SurfGlacierDEM, 'Z', Xbed, Ybed)
    ZSurfDebris = Interpolate_field(Df_SurfDebrisDEM, 'Z', Xbed, Ybed)
    Hdebris_notclipped = ZSurfDebris - ZSurfGlacier
    ### If H debris is less than 1.5m -> Zero
    Hdebris = np.where(Hdebris_notclipped < 1.5, 0, Hdebris_notclipped)
    Hdebris_smoothed = gaussian_filter(Hdebris, sigma=6)

    taub = Interpolate_field(Df_wholebed, 'taub', Xbed, Ybed)
    N = Interpolate_field(Df_wholebed, 'N', Xbed, Ybed)
    ub = Interpolate_field(Df_wholebed, 'ub', Xbed, Ybed)
    Sigma_nn = Interpolate_field(Df_wholebed, 'Sigma_nn', Xbed, Ybed)
    theta = Interpolate_field(Df_wholebed, 'theta', Xbed, Ybed)
    GM = Interpolate_field(Df_wholebed, 'GM', Xbed, Ybed)
    #############################################
    ###     PLOT MAPS OF VARIABLE OF INTEREST ###
    #############################################
    ###~~~~~~~FIG 1 : Theta ~~~~~~~~~~~~~~~
    ###Prepare the subplot
    if VariableToPlot == 'Cavitation':
        fig1, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
        ax.set_xlabel(r'X [km]', fontsize=18)
        ax.set_ylabel(r'Y [km]', fontsize=18)
        ax.tick_params(labelsize=16)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        # ax.set_title('Normal Stress [MPa]', fontsize=21, weight='bold')
        # shading
        clevs = np.arange(0.0, 1.01, 0.1)  ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'cividis'
        cmap = cmx.get_cmap(Colormap)
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, 1 - theta, levels=clevs, cmap=cmap, extend='both', alpha=0.7,antialiased=True)
        # CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, theta, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax, alpha=0.6, zorder=1)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed, levels=DebrisContour_Levels,colors='k', linewidths=DebrisContour_Thick, linestyles='-')
        ###Put labels with position determined manually
        ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d", manual=DebrisContour_LabelPositions, fontsize=7.5)
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ### At floating areas (GM <0), make a contour
        CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.1], colors='black', linewidths=1.2, linestyles='-')
        CS1_hatch = ax.contourf(Xbed / 1000, Ybed / 1000, GM, levels=[-1, -0.1], colors='none', hatches=['////'], alpha=0)
        # ###Plot stars where we plot variable versus time
        if Stars:
            for xstar, ystar, colstar in zip(xplot, yplot, colplot):
                ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Marker_Size,
                        markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1, CS1_debriscontour, CS1_GLcontour, CS1_hatch]:
            c.set_clip_path(patch)

        cbar = fig1.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
        cbar.set_ticks(np.arange(0, 1.01, 0.2))
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label(r'Cavitation', fontsize=28)

    ###~~~~~~~FIG 2 : taub/N ~~~~~~~~~~~~~~~
    if VariableToPlot == 'Frottement':
        ### Make the division by being careful of places where N might be zero or NaN
        taub_over_N_raw = np.divide(taub, N, out=np.full_like(taub, np.nan), where=(N != 0))
        # Mask NaNs and infs for safe plotting
        taub_over_N = np.ma.masked_invalid(taub_over_N_raw)
        ###Prepare the subplot
        fig2, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
        ax.set_xlabel(r'X [km]', fontsize=18)
        ax.set_ylabel(r'Y [km]', fontsize=18)
        ax.tick_params(labelsize=16)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        # ax.set_title('Normal Stress [MPa]', fontsize=21, weight='bold')
        # shading
        clevs = np.arange(0.0, 0.501, 0.05)  ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'YlGnBu_r'
        cmap = cmx.get_cmap(Colormap)
        #colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed/1000, Ybed/1000, taub_over_N, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
        #CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, theta, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax, alpha=0.6, zorder=1)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed, levels=DebrisContour_Levels,colors='k', linewidths=DebrisContour_Thick, linestyles='-')
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ###Put labels with position determined manually
        ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d",manual=DebrisContour_LabelPositions, fontsize=7.5)
        ### At floating areas (GM <0), make a contour
        CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.1], colors='black', linewidths=1.2, linestyles='-')
        # ###Plot stars where we plot variable versus time
        if Stars:
            for xstar, ystar, colstar in zip(xplot, yplot, colplot):
                ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Marker_Size,
                        markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1, CS1_debriscontour, CS1_GLcontour]:
            c.set_clip_path(patch)

        cbar = fig2.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
        cbar.set_ticks(np.arange(0, 0.501, 0.1))
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label(r'Frottement normalisé', fontsize=28)

    ###~~~~~~~FIG 3 : Sigma_nn ~~~~~~~~~~~~~~~
    if VariableToPlot == 'Sigma_nn':
        ###Prepare the subplot
        fig3, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
        ax.set_xlabel(r'X [km]', fontsize=18)
        ax.set_ylabel(r'Y [km]', fontsize=18)
        ax.tick_params(labelsize=16)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        # shading
        clevs = np.arange(-2.0, 2.01, 0.05)  ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'RdBu'
        cmap = cmx.get_cmap(Colormap)
        # colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        # CS1 = ax.contourf(Xbed/1000, Ybed/1000, Sigma_nn, levels=clevs, cmap=cmap, extend='both', alpha=0.7, antialiased=True)
        CS1 = ax.pcolormesh(Xbed / 1000, Ybed / 1000, Sigma_nn, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax, alpha=0.6,zorder=1)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed, levels=DebrisContour_Levels,colors='k', linewidths=DebrisContour_Thick, linestyles='-')
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ###Put labels with position determined manually
        ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d",manual=DebrisContour_LabelPositions, fontsize=7.5)
        ### At floating areas (GM <0), make a contour
        CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.1], colors='black', linewidths=1.2, linestyles='-')
        CS1_hatch = ax.contourf(Xbed / 1000, Ybed / 1000, GM, levels=[-1, -0.1], colors='none', hatches=['////'], alpha=0)
        ###Plot stars where we plot variable versus time
        if Stars:
            for xstar, ystar, colstar in zip(xplot, yplot, colplot):
                ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Marker_Size,
                        markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1, CS1_debriscontour, CS1_GLcontour,CS1_hatch]:
            c.set_clip_path(patch)
        cbar = fig3.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
        cbar.set_ticks(np.arange(-2, 2.01, 0.5))
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label(r'Contrainte Normale [MPa]', fontsize=28)

    ###~~~~~~~FIG 4 : ub ~~~~~~~~~~~~~~~
    if VariableToPlot == 'Velocity':
        ###Prepare the subplot
        fig4, ax = plt.subplots(figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
        ax.set_xlabel(r'X [km]', fontsize=18)
        ax.set_ylabel(r'Y [km]', fontsize=18)
        ax.tick_params(labelsize=16)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        # shading
        clevs = np.arange(0.0, velocitymax + 0.1, 2.5)  ## cbar for shading
        vmin, vmax = clevs[0], clevs[-1]
        Colormap = 'viridis'
        cmap = cmx.get_cmap(Colormap)
        # colorbar
        # levs_ticks = [1, 10, 100, 1000, 10000]
        # --- Plot orthophoto as background
        ax.imshow(ortho_crop, extent=extent_km, origin='upper', cmap='gray', zorder=0)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(Xbed / 1000, Ybed / 1000, ub / 365, levels=clevs, cmap=cmap, extend='both', alpha=0.7,antialiased=True)
        ### Plot contours of debris thickness
        CS1_debriscontour = ax.contour(Xbed / 1000, Ybed / 1000, Hdebris_smoothed, levels=DebrisContour_Levels,colors='white', linewidths=DebrisContour_Thick, linestyles='-')
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ###Put labels with position determined manually
        ax.clabel(CS1_debriscontour, DebrisContour_Levels, inline=True, inline_spacing=2, fmt="%d",manual=DebrisContour_LabelPositions, fontsize=7.5)
        ### At floating areas (GM <0), make a contour
        CS1_GLcontour = ax.contour(Xbed / 1000, Ybed / 1000, GM, levels=[-0.1], colors='black', linewidths=1.2, linestyles='-')
        # Hatching the area where GM < -0.1
        # We set vmin, vmax so that contourf only fills below -0.1
        CS1_hatch = ax.contourf(Xbed / 1000, Ybed / 1000, GM,levels=[-1, -0.1],  colors='none',  hatches=['////'],  alpha=0 )
        ###Plot stars where we plot variable versus time
        if Stars:
            for xstar, ystar, colstar in zip(xplot, yplot, colplot):
                ax.plot(xstar / 1000, ystar / 1000, color=colstar, marker='*', markersize=Marker_Size,
                        markeredgecolor='black', markeredgewidth=1.2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1, CS1_debriscontour, CS1_GLcontour, CS1_hatch]:
            c.set_clip_path(patch)
        cbar = fig4.colorbar(CS1, ax=ax, shrink=1, pad=0.05)
        cbar.set_ticks(np.arange(0, velocitymax + 0.1, 5))
        cbar.ax.tick_params(labelsize=18)
        cbar.set_label(r'Vitesse basale [m/jour]', fontsize=28)




    plt.show()

    # ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    # ####      COMMON TO THE TWO SUBPLOTS        ####
    # ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    # Add a common colorbar to all subplots below the figure
    # cbar_ax = fig1.add_axes([0.15, 0.082, 0.7, 0.027])  # [left, bottom, width, height]
    # fig1.colorbar(CS1,  cax=cbar_ax, orientation='horizontal', label=r'$\tau_\mathrm{b}$/N')
    # fig1.colorbar(CS1, ticks=levs_ticks, cax=cbar_ax, orientation='horizontal', label=r'Maxwell time [h]')
    # cbar_ax.set_xticklabels([r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
    # ###Show map
    # plt.show()




