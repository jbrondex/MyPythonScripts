"""
@author: jbrondex

Description:
------------
This file aims at creating a geo (gmsh input file) from the 2D (x,y) DEM file bed corresponding to 7 transects extracted from 3D DEMs of bed of Grande Motte Glacier (Tignes)
using QGIS. Here the purpose is simply to create a mesh with bed DEM as top surface and random flat elevation as bottom surface to have
a vtu that enables to easily visualize the bed on paraview

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
el_size = 0.1
el_sizec = 0.1
# Spline or line
spline = True ###It seems to be working all the time here so True by default
#### Transect and corresponding tunnel floor altitude (1% slope)
Transect_Names = ['Transect1','Transect2','Transect3','Transect4','Transect5','Transect6']
Floor_altitude= 2770

#########################################################
########## START LOOP ON TRANSECTS ######################
#########################################################
for Transect in Transect_Names:
    ###Get the bed DEM of the domain contour
    file_name_bed = './DEMs/DEM_Bed_{}.dat'.format(Transect)
    # Load DEM of domain
    Col_Names = ['x', 'y']
    DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
    Nptbed = len(DEM_bed.index)
    print('Transect n°', Transect)
    file_name_output = 'GrandeMotte_BedVisu_{}.geo'.format(Transect)
    geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
    geo.write('// This a a geo file created using the python script Makegeo_GrandeMotteVisuBed_Transect.py // \n')
    geo.write('Mesh.Algorithm=5; \n')
    geo.write('// To control the element size, one can directly modify the lc values in the geo file // \n')
    geo.write('lc = {} ; \n'.format(el_size))
    geo.write('lcc = {} ; \n'.format(el_sizec))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # write the points coordinates (x,y,0,lc) FOR THE DOMAIN CONTOUR
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###DEM top and bed have to be treated separately to have their own BCs for the domain contour
    np=0
    DomainTopLeftCornerPoint_Index = np + 1
    ###...continu with bottom nodes
    for j in range(0,Nptbed):
        np=np+1
        geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_bed['x'][j],DEM_bed['y'][j])+r'}'+'; \n')
    ###Store the index of Domain bottom left corner
    DomainTopRightCornerPoint_Index = np
    ###... continue with bottom nodes
    np = np + 1
    DomainBottomRightCornerPoint_Index = np
    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(DEM_bed['x'].max(), Floor_altitude) + r'}' + '; \n')
    np = np + 1
    DomainBottomLeftCornerPoint_Index = np
    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(DEM_bed['x'].min(), Floor_altitude) + r'}' + '; \n')
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###Now build the lines from the points
    # if spline
    if spline:
        ###The spline corresponding to the top surface of the domain
        geo.write('Spline(1) = {')
        for j in range(1,DomainTopRightCornerPoint_Index):
            geo.write('{},'.format(j))
        geo.write('{}'.format(DomainTopRightCornerPoint_Index)+'}; \n')
        ###The spline corresponding to right side of the domain
        geo.write('Spline(2) = {')
        geo.write('{},'.format(DomainTopRightCornerPoint_Index))
        geo.write('{}'.format(DomainBottomRightCornerPoint_Index)+'}; \n')
        ###The spline corresponding to the bed of the domain
        geo.write('Spline(3) = {')
        geo.write('{},'.format(DomainBottomRightCornerPoint_Index))
        geo.write('{}'.format(DomainBottomLeftCornerPoint_Index) + '}; \n')
        ###The spline corresponding to left side of the domain
        geo.write('Spline(4) = {')
        geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
        geo.write('1}; \n')
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ###Now build the line loops
        geo.write('Line Loop(5) = {1, 2, 3, 4}; \n') ###domain contour
        geo.write('Plane Surface(6) = {5}; \n')
        geo.write('Physical Line(7) = {1}; \n') ##BC corresponding to bed DEM (I used it to export node coordinates)
        geo.write('Physical Line(8) = {2, 3, 4}; \n') ##Other BCs
        geo.write('Physical Surface(9) = {6}; \n')
    geo.close()

