"""
@author: jbrondex

Description:
------------
This file aims at creating a geo (gmsh input file) from two 2D (x,y) DEM files (top and bed) corresponding to 7 transects extracted from 3D DEMs of surface and bed of Grande Motte Glacier (Tignes)
using QGIS. For each transect we add a tunnel footprint considering 3 possible shapes: rectangle, rectangle + half-circle, ovoide.

"""

# -*- coding: utf-8 -*-
# Create a geo (gmsh input file) file from two contour files
# the contour files contain the (x,y) coordinates of the ordered
# points defining the contour of the glacier and the contour of the cavity
#
import numpy as np
import pandas as pd
from pathlib import Path

### Where to find the DEM files and store the .geo file:
pathroot_mycode = Path('/home/brondexj/Documents/GrandeMotteTignes/Maillages/.')

###################################################
######## Parameters of the various meshes #########
###################################################
# edge size of the elements
el_size = 0.6
el_sizec = 0.1
# Spline or line
spline = True ###It seems to be working all the time here so True by default
#### Transect and corresponding tunnel floor altitude (1% slope)
Transect_Names = ['Transect1','Transect2','Transect3','Transect4','Transect5','Transect6']
Floor_altitudes=[2803.5,2803,2802.5,2802,2801.5,2801.12]
#### Possible shapes of the tunnel
Shapes = ['Rectangle', 'Circle', 'Ovoide']
#### Width of tunnel (for Rectangle and Half Circle case only) and corresponding name (in cm) for output file name
Widths = [1,1.5,2]
Width_Names=['100','150','200']

###Build the ovoid DEM from parameterized curve method for the ovoid shape (based on RTM ppt)
t = np.linspace(-5, 5, 301)
x = 3.2 * t / (1 + t ** 2) ** 2
y = 3.2 / (1 + t ** 2) ** 2
df = pd.DataFrame({'xo': x, 'yo': y})
#### Remove points that are below 2.5m of maximum height
df = df[df['yo'] > np.max(df['yo']) - 2.5]
###Add a last point corresponding to the center of the floor as the last line
df.loc[len(df.index)] = [0.0, np.max(df['yo']) - 2.5]  # adding a row
df = df.reset_index(drop=True)
Floor_Ovoid = np.min(df['yo'])
#########################################################
########## START LOOP ON TRANSECTS ######################
#########################################################
for (Transect, Floor_altitude) in zip(Transect_Names,Floor_altitudes):
    ###Get the bed and top DEM of the domain contour
    file_name_top = './DEMs/DEM_Top_{}.dat'.format(Transect)
    file_name_bed = './DEMs/DEM_Bed_{}.dat'.format(Transect)
    # Load DEM of domain
    Col_Names = ['x', 'y']
    DEM_top = pd.read_csv(pathroot_mycode.joinpath(file_name_top), names=Col_Names, delim_whitespace=True)
    Npttop = len(DEM_top.index)
    DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
    ###Reverse DEM_bed dataframe so that nodes are consecutive
    DEM_bed = DEM_bed[::-1]
    DEM_bed= DEM_bed.reset_index(drop=True)
    Nptbed = len(DEM_bed.index)
    ###Merge both DEM to get full contour of domain
    cont=[DEM_top, DEM_bed]
    MainContour = pd.concat(cont)
    MainContour=MainContour.reset_index(drop=True)

    ###################################################
    ########## START LOOP ON SHAPES ###################
    ###################################################
    for Shape in Shapes:
        print('Transect n°', Transect, '   Shape:', Shape)
        ###################################################
        ########## START LOOP ON WIDTHS ###################
        ###################################################
        for k,(Width,Width_Name) in enumerate(zip(Widths, Width_Names)):
            if Shape == 'Ovoide':
                if k==0:
                    ###For the ovoid case, the tunnel must be centered based on floor altitude
                    DEM_Ovoid = pd.DataFrame({'xt': df['xo'] + 25, 'yt': Floor_altitude - Floor_Ovoid + df['yo']})
                else: ###for the ovoid tunnel, we do not try several widths so only go through the loop once
                    break
            #######################################
            ###NOW START WRITING THE .geo FILE ####
            #######################################
            # Open the output file and write some general info
            if Shape == 'Ovoide':
                ####Set the name of the output file
                file_name_output = 'GrandeMotte_{}_{}.geo'.format(Transect, Shape)
            else:
                file_name_output='GrandeMotte_{}_{}_W{}cm.geo'.format(Transect,Shape,Width_Name)
            geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
            geo.write('// This a a geo file created using the python script Makegeo_GrandeMotteTunnel_Transect.py // \n')
            geo.write('Mesh.Algorithm=5; \n')
            geo.write('// To control the element size, one can directly modify the lc values in the geo file // \n')
            geo.write('lc = {} ; \n'.format(el_size))
            geo.write('lcc = {} ; \n'.format(el_sizec))
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # write the points coordinates (x,y,0,lc) FOR THE DOMAIN CONTOUR
            ###DEM top and bed have to be treated separately to have their own BCs for the domain contour
            np=0
            ###start by top nodes ...
            for j in range(0,Npttop):
                np=np+1
                geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_top['x'][j],DEM_top['y'][j])+r'}'+'; \n')
            ###...continu with bottom nodes
            for j in range(0,Nptbed):
                np=np+1
                geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_bed['x'][j],DEM_bed['y'][j])+r'}'+'; \n')
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ###...now write the coordinates (x,y,0,lc) FOR THE TUNNEL CONTOUR
            ####For the rectangle and half circle case, only a few points are required (it is important to turn counter-clokwise for the half-circle case)
            if Shape == 'Rectangle' or Shape == 'Circle':
                ##Bottom left corner
                geo.write('Point({}) = '.format(np+1) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude) + r'}' + '; \n')
                ###Bottom right corner
                geo.write('Point({}) = '.format(np+2) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude) + r'}' + '; \n')
                ###Top corners
                if Shape == 'Rectangle':
                    ### Top right
                    geo.write('Point({}) = '.format(np+3) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5) + r'}' + '; \n')
                    ### Top Left
                    geo.write('Point({}) = '.format(np+4) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5) + r'}' + '; \n')
                elif Shape == 'Circle':
                    ### Top right
                    geo.write('Point({}) = '.format(np+3) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
                    ### Top Left
                    geo.write('Point({}) = '.format(np+4) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
                ###If circle, requires circle center as well
                if Shape == 'Circle':
                    geo.write('Point({}) = '.format(np+5) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
            ####For the ovoid shape more points are required
            elif Shape == 'Ovoide':
                ####Get points coordinates from the DEM dataframe
                for j in range(len(DEM_Ovoid)):
                    np = np + 1
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(DEM_Ovoid['xt'][j], DEM_Ovoid['yt'][j]) + r'}' + '; \n')
            ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ###Now build the lines from the points
            # if spline
            if spline:
                ###The spline corresponding to the top surface of the domain
                geo.write('Spline(1) = {')
                for j in range(0,Npttop-1):
                    geo.write('{},'.format(j+1))
                geo.write('{}'.format(Npttop)+'}; \n')
                ###The spline corresponding to right side of the domain
                geo.write('Spline(2) = {')
                geo.write('{},'.format(Npttop))
                geo.write('{}'.format(Npttop+1)+'}; \n')
                ###The spline corresponding to the bed of the domain
                geo.write('Spline(3) = {')
                for j in range(0,Nptbed-1):
                    geo.write('{},'.format(Npttop+j+1))
                geo.write('{}'.format(Npttop + Nptbed) + '}; \n')
                ###The spline corresponding to right side of the domain
                geo.write('Spline(4) = {')
                geo.write('{},'.format(Npttop+Nptbed))
                geo.write('1}; \n')
                ###Now the spline for the tunnel (important to turn counter-clockwise for the half circle case)
                if Shape == 'Rectangle' or Shape == 'Circle':
                    ###Bottom line
                    geo.write('Spline(5) = {')
                    geo.write('{},'.format(Npttop+Nptbed + 1))
                    geo.write('{}'.format(Npttop + Nptbed + 2))
                    geo.write('}; \n')
                    ###Right line
                    geo.write('Spline(6) = {')
                    geo.write('{},'.format(Npttop+Nptbed + 2))
                    geo.write('{}'.format(Npttop + Nptbed + 3))
                    geo.write('}; \n')
                    ###Top line (or half circle)
                    if Shape == 'Rectangle':
                        geo.write('Spline(7) = {')
                        geo.write('{},'.format(Npttop + Nptbed + 3))
                        geo.write('{}'.format(Npttop + Nptbed + 4))
                        geo.write('}; \n')
                    elif Shape == 'Circle':
                        geo.write('Circle(7) = {')
                        geo.write('{},'.format(Npttop + Nptbed + 3))
                        geo.write('{},'.format(Npttop + Nptbed + 5))
                        geo.write('{}'.format(Npttop + Nptbed + 4))
                        geo.write('}; \n')
                    geo.write('Spline(8) = {')
                    geo.write('{},'.format(Npttop + Nptbed + 4))
                    geo.write('{}'.format(Npttop + Nptbed + 1))
                    geo.write('}; \n')
                elif Shape == 'Ovoide':
                    ###Ovoid walls without floor
                    geo.write('Spline(5) = {')
                    for j in range(len(DEM_Ovoid)-2):
                        geo.write('{},'.format(Npttop+Nptbed + j + 1))
                    geo.write('{}'.format(Npttop+Nptbed + len(DEM_Ovoid)-1) + '}; \n')
                    ###Ovoid floor
                    geo.write('Spline(6) = {')
                    geo.write('{},'.format(Npttop + Nptbed + len(DEM_Ovoid) - 1))
                    geo.write('{},'.format(Npttop + Nptbed + len(DEM_Ovoid)))
                    geo.write('{}'.format(Npttop + Nptbed + 1))
                    geo.write('}; \n')
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ###Now build the line loops
                geo.write('Line Loop(9) = {1, 2, 3, 4}; \n') ###domain contour
                if Shape == 'Rectangle' or Shape == 'Circle': ###tunnel contour
                    geo.write('Line Loop(10) = {5, 6, 7, 8}; \n')
                elif Shape == 'Ovoide':
                    geo.write('Line Loop(10) = {5, 6}; \n')
                geo.write('Plane Surface(11) = {9,10}; \n')
                geo.write('Physical Line(12) = {1}; \n') ##BC top
                geo.write('Physical Line(13) = {2}; \n') ##BC right side
                geo.write('Physical Line(14) = {3}; \n') ##BC bottom
                geo.write('Physical Line(15) = {4}; \n') ##BC left side
                geo.write('Physical Line(16) = {5, 6, 7, 8}; \n') ##BC tunnel (Free surface)
                geo.write('Physical Surface(17) = {11}; \n')
            geo.close()

