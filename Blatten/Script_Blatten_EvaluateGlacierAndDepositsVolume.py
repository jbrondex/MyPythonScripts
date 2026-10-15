"""
@author: jbrondex

Description:
------------
This file aims at evaluating initial cavity volume

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
from scipy.integrate import simpson

from matplotlib import ticker
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


if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###For which Year do we want to make the calculation
    Year = 2025
    ###Open Data corresponding to DEM bed and bottom surf
    Col_Names = ['X', 'Y', 'Z']
    # Pathroot = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/1_Mesh/MESH2025/.')
    Pathroot = Path('/home/brondexj/Documents/Expertises/Blatten/DEMs/.')
    Filename_zbRock = 'DEM_surf_2025_Reconstructed_FromAdrien_MJacquemart_Smoothed.dat' ##'DEM_surf_Corrected.dat'
    Filename_zsRock = 'DEMSurfPre.dat'
    Filename_zbed = 'DEM_bed_Corrected.dat'
    ###load files as dataframes
    Df_zbRock = pd.read_csv(Pathroot.joinpath(Filename_zbRock), names=Col_Names, delim_whitespace=True)
    Df_zsRock = pd.read_csv(Pathroot.joinpath(Filename_zsRock), names=Col_Names, delim_whitespace=True)
    Df_zbed = pd.read_csv(Pathroot.joinpath(Filename_zbed), names=Col_Names, delim_whitespace=True)
    ###For the DEM of the glacier surface, name depends on considered year
    #Pathrootzs = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup_PROD_FINAL_Mai2026/1_Mesh/MESH{}/.'.format(Year))
    Pathrootzs = Path('/home/brondexj/Documents/Expertises/Blatten/DEMs/.')
    if Year == 2025:
        Filename_zsGlacier = 'DEM_surf_2025_Reconstructed_FromAdrien_MJacquemart_Smoothed.dat' ##'DEM_surf_Corrected.dat'
    else:
        Filename_zsGlacier = 'DEMSurf{}.dat'.format(Year)
    Df_zsGlacier = pd.read_csv(Pathrootzs.joinpath(Filename_zsGlacier), names=Col_Names, delim_whitespace=True)
    ###load contour of glacier (to adapt depending on Year)
    Filename_GlacierContour = 'Contour_Birch_{}_Julien.dat'.format(Year)
    ###Print some info
    print('Opening file:', Pathrootzs.joinpath(Filename_GlacierContour))
    ###Load glacier contour
    Df_GlacierContour = pd.read_csv(Pathrootzs.joinpath(Filename_GlacierContour), names=['X', 'Y'], sep=r'\s+')
    # Close the contour: Copy the first row and append it to the end
    Df_GlacierContour = pd.concat([Df_GlacierContour, Df_GlacierContour.iloc[[0]]], ignore_index=True)
    ### Get xmin, xmax, ymin, ymax of contour
    xmin = Df_GlacierContour['X'].min()
    xmax = Df_GlacierContour['X'].max()
    ymin = Df_GlacierContour['Y'].min()
    ymax = Df_GlacierContour['Y'].max()
    ### Create a regular grid common to all DEMs
    dx = 0.25
    dy = 0.25
    xgrid= np.arange(xmin, xmax, dx)
    ygrid= np.arange(ymin, ymax, dy)
    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    Zbed = Interpolate_field(Df_zbed, 'Z', Xgrid, Ygrid)
    ZbRock = Interpolate_field(Df_zbRock, 'Z', Xgrid, Ygrid)
    ZsRock = Interpolate_field(Df_zsRock, 'Z', Xgrid, Ygrid)
    ZsGlacier = Interpolate_field(Df_zsGlacier, 'Z', Xgrid, Ygrid)
    ###Force all z of all DEMs to zero for DEM points that are out of the glacier contour
    # 1. Build a Path object from the contour
    poly_path = mpltPath(np.column_stack((Df_GlacierContour['X'], Df_GlacierContour['Y'])))
    # 2. Flatten your grid points (Xgrid, Ygrid) -> (N, 2)
    points = np.column_stack((Xgrid.ravel(), Ygrid.ravel()))
    # 3. Check which grid points are inside the polygon
    inside_mask = poly_path.contains_points(points, radius=-1e-6)
    # 4. Reshape mask back to 2D grid shape
    inside_mask = inside_mask.reshape(Xgrid.shape)
    # 5. Apply mask to all your Z grids
    Zbed_masked = np.where(inside_mask, Zbed, 0)
    ZbRock_masked = np.where(inside_mask, ZbRock, 0)
    ZsRock_masked = np.where(inside_mask, ZsRock, 0)
    ZsGlacier_masked = np.where(inside_mask, ZsGlacier, 0)

    #### COMPUTE VOLUME OF DEBRIS AND GLACIER ####
    ###Calculate debris thickness
    Hrock = ZsRock_masked - ZbRock_masked
    # Force negative thickness to zero
    Hrock[Hrock < 0] = 0
    # total volume
    VolumeRock = np.sum(Hrock) * dx * dy  # in cubic meters
    print(f"Max thickness Rock: {Hrock.max():.2f} m")
    print(f"Mean thickness Rock (non-zero): {Hrock[Hrock > 0].mean():.2f} m")
    print(f"Number of non-zero cells Rock: {(Hrock > 0).sum()}")
    print(f"Total rock volume: {VolumeRock/1e6:.2f} millions of m³")


    ###Same thing for glacier
    Hglacier = ZsGlacier_masked - Zbed_masked
    # Force negative thickness to zero
    Hglacier[Hglacier< 0] = 0
    # total volume
    VolumeGlacier = np.sum(Hglacier) * dx * dy  # in cubic meters
    print(f"Max thickness Glacier: {Hglacier.max():.2f} m")
    print(f"Mean thickness Glacier (non-zero): {Hglacier[Hglacier > 0].mean():.2f} m")
    print(f"Number of non-zero cells Glacier: {(Hglacier > 0).sum()}")
    print(f"Total glacier volume: {VolumeGlacier/1e6:.2f} millions of m³")

