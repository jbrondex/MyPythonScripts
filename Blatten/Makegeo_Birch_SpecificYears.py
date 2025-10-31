"""
@author: jbrondex

Description:
------------
This file aims at creating a geo (gmsh input file) for the Birch glacier for specific years

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
pathroot_mycode = Path('/home/brondexj/BETTIK/BlattenEventAdri/MyRateAndStateSetup/1_Mesh/.')

###################################################
######## Parameters of the various meshes #########
###################################################
# edge size of the elements
el_size = 10.0
#el_sizec = 0.1
# Spline or line
spline = True ###It seems to be working all the time here so True by default
#### Years for which we want a mesh and corresponding last node of side BC (following node correspond to front boundary)
Years= ['1993','2011']
LastNodeSide=[38,35]

#########################################################
########## START LOOP ON YEARS ######################
#########################################################
for (Year, LastNode) in zip(Years,LastNodeSide):
    ###Get the domain contour
    file_name_cont = './MESH{}/Contour_Birch_{}.dat'.format(Year,Year)
    # Load DEM of domain
    Col_Names = ['x', 'y']
    Contour = pd.read_csv(pathroot_mycode.joinpath(file_name_cont), names=Col_Names, sep=r'\s+')
    NptTot = len(Contour.index)
    ### Declare Name of Output
    file_name_output = './MESH{}/Contour_Birch_{}.geo'.format(Year,Year)

    #######################################
    ###NOW START WRITING THE .geo FILE ####
    #######################################
    geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
    geo.write('// This a a geo file created using the python script Makegeo_Birch_SpecificYears.py // \n')
    geo.write('Mesh.Algorithm=5; \n')
    geo.write('// To control the element size, one can directly modify the lc values in the geo file // \n')
    geo.write('lc = {} ; \n'.format(el_size))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # write the points coordinates (x,y,0,lc) FOR THE DOMAIN CONTOUR
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###Nodes of side and front have to be treated separately to have their own BCs for the domain contour
    np=0
    ###Points
    for j in range(0,NptTot):
        np=np+1
        geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(Contour['x'][j],Contour['y'][j])+r'}'+'; \n')
    ###Now build the lines from the points
    # if spline
    if spline:
        ###The spline corresponding to the side
        geo.write('Spline(1) = {')
        for k in range(1,LastNode):
            geo.write('{},'.format(k))
        geo.write('{}'.format(LastNode)+'}; \n')
        ###The spline corresponding to the front
        geo.write('Spline(2) = {')
        for k in range(LastNode,NptTot+1):
            geo.write('{},'.format(k))
        geo.write('1}; \n')
        ###Now build the line loops
        geo.write('Line Loop(3) = {1,2}; \n') ###domain contour
        geo.write('Plane Surface(4) = {3}; \n')
        geo.write('Physical Line(5) = {1}; \n') ##BC side
        geo.write('Physical Line(6) = {2}; \n') ##BC front
        geo.write('Physical Surface(7) = {4}; \n')
    geo.close()

