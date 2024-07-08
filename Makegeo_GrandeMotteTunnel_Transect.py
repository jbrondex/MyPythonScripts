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

# Test these options
    # edge size of the elements
el_size = 0.8
el_sizec = 0.1

    # Spline or line 
spline = True

    # Load DEM
Col_Names = ['x', 'y']
DEM_top = pd.read_csv('./DEM_Top_Transect1.dat', names=Col_Names, delim_whitespace=True)
Npttop = len(DEM_top.index)
DEM_bed = pd.read_csv('./DEM_Bed_Transect1.dat', names=Col_Names, delim_whitespace=True)
###Reverse DEM_bed dataframe so that nodes are consecutive
DEM_bed = DEM_bed[::-1]
DEM_bed= DEM_bed.reset_index(drop=True)
Nptbed = len(DEM_bed.index)
###Merge both DEM to get full contour
cont=[DEM_top, DEM_bed]
MainContour = pd.concat(cont)
MainContour=MainContour.reset_index(drop=True)


# Open the output file
geo = open('GrandeMotte_Transect1_Rectangle.geo', 'w')
geo.write('// This a a geo file created using the python script Makegeo_2.py // \n')
geo.write('Mesh.Algorithm=5; \n')
geo.write('// To controle the element size, one can directly modify the lc values in the geo file // \n')
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
    geo.write('Point({0}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_bed['x'][j],DEM_bed['y'][j])+r'}'+'; \n')

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
    
    geo.write('Line Loop(5) = {1, 2, 3, 4}; \n')
#    geo.write('Line Loop(4) = {2}; \n')
#    geo.write('Plane Surface(6) = {3,4}; \n')
    geo.write('Plane Surface(6) = {5}; \n')
    geo.write('Physical Line(7) = {1}; \n')
    geo.write('Physical Line(8) = {2}; \n')
    geo.write('Physical Line(9) = {3}; \n')
    geo.write('Physical Line(10) = {4}; \n')
    geo.write('Physical Surface(11) = {6}; \n')
    
    
# else it is lines, as a spline might not work in all case
else:
    nl=0
    for j in range(0,Npt-1):
        nl=nl+1
        geo.write('Line({0}) = '.format(nl)+r'{'+'{0},{1}'.format(j+1,j+2)+r'}'+'; \n')
    geo.write('Line({0}) = '.format(nl+1)+r'{'+'{0},{1}'.format(j+2,1)+r'}'+'; \n')
    
    nl = nl+1
    for j in range(0,Nptc-1):
        nl=nl+1
        geo.write('Line({0}) = '.format(nl)+r'{'+'{0},{1}'.format(Npt+j+1,Npt+j+2)+r'}'+'; \n')
    geo.write('Line({0}) = '.format(nl+1)+r'{'+'{0},{1}'.format(Npt+j+2,Npt+1)+r'}'+'; \n')
    
    geo.write('Compound Line({0}) = '.format(nl+2)+r'{')
    for j in range(0,Npt-1):
        geo.write('{0}, '.format(j+1))
    geo.write('{0}'.format(j+2)+'}; \n')
    
    geo.write('Compound Line({0}) = '.format(nl+3)+r'{')
    for j in range(0,Nptc-1):
        geo.write('{0}, '.format(Npt+j+1))
    geo.write('{0}'.format(Npt+j+2)+'}; \n')
    
    geo.write('Line Loop({0}) = '.format(nl+4)+r'{'+'{0}'.format(nl+2)+r'};'+' \n')
    geo.write('Line Loop({0}) = '.format(nl+5)+r'{'+'{0}'.format(nl+3)+r'};'+' \n')
    
    geo.write('Plane Surface({0}) = '.format(nl+6)+r'{'+'{0},{1}'.format(nl+4,nl+5)+r'};'+' \n')
    geo.write('Plane Surface({0}) = '.format(nl+7)+r'{'+'{0}'.format(nl+5)+r'};'+' \n')
    
    geo.write('Physical Line({0}) = '.format(nl+8)+r'{'+'{0}'.format(nl+2)+r'};'+' \n')
    geo.write('Physical Surface({0}) = '.format(nl+9)+r'{'+'{0},{1}'.format(nl+6,nl+7)+r'};'+' \n')

geo.close()
