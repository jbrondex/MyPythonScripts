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
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

###Parameters of the various meshes
### Where to find the DEM files and store the .geo file:
pathroot_mycode = Path('/home/brondexj/Documents/GrandeMotteTignes/Maillages/.')
# Test these options
    # edge size of the elements
el_size = 0.6
el_sizec = 0.1
    # Spline or line 
spline = True
    #### Transect and corresponding tunnel floor altitude (1% slope)
Transect_Names = ['Transect1','Transect2','Transect3','Transect4','Transect5','Transect6']
Transect=Transect_Names[0]
Floor_altitudes=[2803.5,2803,2802.5,2802,2801.5,2801.12]
Floor_altitude=Floor_altitudes[0]
    #### Possible shapes of the tunnel
Shapes = ['Rectangle', 'Circle', 'Ovoide']
Shape= Shapes[0]
    #### Width of tunnel
Widths = [1,1.5,2]
Width = Widths[2]
Width_Names=['100','150','200']
Width_Name = Width_Names[2]

file_name_top = 'DEM_Top_{}.dat'.format(Transect)
file_name_bed = 'DEM_Bed_{}.dat'.format(Transect)
    # Load DEM
Col_Names = ['x', 'y']
DEM_top = pd.read_csv(pathroot_mycode.joinpath(file_name_top), names=Col_Names, delim_whitespace=True)
Npttop = len(DEM_top.index)
DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
###Reverse DEM_bed dataframe so that nodes are consecutive
DEM_bed = DEM_bed[::-1]
DEM_bed= DEM_bed.reset_index(drop=True)
Nptbed = len(DEM_bed.index)
###Merge both DEM to get full contour
cont=[DEM_top, DEM_bed]
MainContour = pd.concat(cont)
MainContour=MainContour.reset_index(drop=True)


# Open the output file
file_name_output='GrandeMotte_{}_{}_W{}cm.geo'.format(Transect,Shape,Width_Name)
geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
geo.write('// This a a geo file created using the python script Makegeo_2.py // \n')
geo.write('Mesh.Algorithm=5; \n')
geo.write('// To control the element size, one can directly modify the lc values in the geo file // \n')
geo.write('lc = {} ; \n'.format(el_size))
geo.write('lcc = {} ; \n'.format(el_sizec))

# write the points coordinates (x,y,0,lc)
###DEM top and bed have to be treated separately to have their own BCs
np=0
###start by top nodes ...
for j in range(0,Npttop):
    np=np+1
    geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_top['x'][j],DEM_top['y'][j])+r'}'+'; \n')
###...continu with bottom nodes
for j in range(0,Nptbed):
    np=np+1
    geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_bed['x'][j],DEM_bed['y'][j])+r'}'+'; \n')
###...now we car about the points corresponding to the tunnel (it is important to turn unclokwise for the half-circle case)
np=np+1
##Bottom left corner
geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude) + r'}' + '; \n')
###Bottom right corner
geo.write('Point({}) = '.format(np+1) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude) + r'}' + '; \n')
###Top corners
if Shape == 'Rectangle':
    ### Top right
    geo.write('Point({}) = '.format(np+2) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5) + r'}' + '; \n')
    ### Top Left
    geo.write('Point({}) = '.format(np+3) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5) + r'}' + '; \n')
elif Shape == 'Circle':
    ### Top right
    geo.write('Point({}) = '.format(np+2) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
    ### Top Left
    geo.write('Point({}) = '.format(np+3) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
###If circle, requires circle center as well
if Shape == 'Circle':
    geo.write('Point({}) = '.format(np+4) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')

# if spline
if spline:
    ###The spline corresponding to the top surface
    geo.write('Spline(1) = {')
    for j in range(0,Npttop-1):
        geo.write('{},'.format(j+1))
    geo.write('{}'.format(Npttop)+'}; \n')
    ###The spline corresponding to right side
    geo.write('Spline(2) = {')
    geo.write('{},'.format(Npttop))
    geo.write('{}'.format(Npttop+1)+'}; \n')
    ###The spline corresponding to the bed
    geo.write('Spline(3) = {')
    for j in range(0,Nptbed-1):
        geo.write('{},'.format(Npttop+j+1))
    geo.write('{}'.format(Npttop + Nptbed) + '}; \n')
    ###The spline corresponding to right side
    geo.write('Spline(4) = {')
    geo.write('{},'.format(Npttop+Nptbed))
    geo.write('1}; \n')
    ###Now the spline for the tunnel (important to turn unclockwise)
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
    
    geo.write('Line Loop(9) = {1, 2, 3, 4}; \n') ###domain contour
    geo.write('Line Loop(10) = {5, 6, 7, 8}; \n') ###tunnel contour
    geo.write('Plane Surface(11) = {9,10}; \n')
    geo.write('Physical Line(12) = {1}; \n') ##BC top
    geo.write('Physical Line(13) = {2}; \n') ##BC right side
    geo.write('Physical Line(14) = {3}; \n') ##BC bottom
    geo.write('Physical Line(15) = {4}; \n') ##BC left side
    geo.write('Physical Line(16) = {5, 6, 7, 8}; \n') ##BC tunnel (Free surface)
    geo.write('Physical Surface(17) = {11}; \n')


geo.close()
