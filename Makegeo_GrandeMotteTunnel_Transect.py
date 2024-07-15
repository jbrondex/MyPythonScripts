"""
@author: jbrondex

Description:
------------
This file aims at creating a geo (gmsh input file) from two 2D (x,y) DEM files (top and bed) corresponding to 7 transects extracted from 3D DEMs of surface and bed of Grande Motte Glacier (Tignes)
using QGIS. For each transect we add a tunnel footprint considering 3 possible shapes: rectangle, rectangle + half-circle, ovoide.
In this script, we can create meshes corresponding to initial shape of tunnel (AfterIncision = False) or to the shape of the tunnel after
incision assuming that the minimum floor altitude of the tunnel is 2792.2 m (-1% slope) which corresponds to the lake surface elevation after draining operations (i.e. if the channel does not reach the bed).
Note that if we assume that water incision as reached the bed on the whole length of the tunnel, corresponding meshes must be created using
the script "Makegeo_GrandeMotteChannel_Transect.py" instead of this one. For the incised tunnel case, there is the possibility, through
the keyword "Incurved_Walls" of having tunnel walls that are not straight but slightly incurved to evaluate wether butressing effects
are able to slow down the tunnel closure.

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
##Do we want tunnel initial shape (False) or tunnel after water incision (True) ?
AfterIncision = False
## Do we want straight walls or incurved walls ?
Incurved_Walls = False
xwallshift = 1.0 ###Controls amount of bending of walls (i.e. the wall middle point is shift by xwallshift relative to straight line)
#### Transect and corresponding tunnel floor altitude (1% slope)
Transect_Names = ['Transect1','Transect2','Transect3','Transect4','Transect5','Transect6']
Floor_altitudes=[2803.5,2803,2802.5,2802,2801.5,2801.12]
BottomTunnel_altitudes=[2792.2,2791.7,2791.2,2790.7,2790.2,2789.82] ### The elevation of tunnel bottom if Tunnel_After_Incision = True
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
for (Transect, Floor_altitude, BottomTunnel_altitude) in zip(Transect_Names,Floor_altitudes,BottomTunnel_altitudes):
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
                if AfterIncision:
                    if Incurved_Walls:
                        file_name_output = 'GrandeMotte_TunnelIncised_IncurvedWalls_{}_{}.geo'.format(Transect, Shape)
                    else:
                        file_name_output = 'GrandeMotte_TunnelIncised_StraightWalls_{}_{}.geo'.format(Transect, Shape)
                else:
                    file_name_output = 'GrandeMotte_TunnelInit_{}_{}.geo'.format(Transect, Shape)
            else:
                if AfterIncision:
                    if Incurved_Walls:
                        file_name_output='GrandeMotte_TunnelIncised_IncurvedWalls_{}_{}_W{}cm.geo'.format(Transect,Shape,Width_Name)
                    else:
                        file_name_output='GrandeMotte_TunnelIncised_StraightWalls_{}_{}_W{}cm.geo'.format(Transect,Shape,Width_Name)
                else:
                        file_name_output='GrandeMotte_TunnelInit_{}_{}_W{}cm.geo'.format(Transect,Shape,Width_Name)

            geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
            geo.write('// This a a geo file created using the python script Makegeo_GrandeMotteTunnel_Transect_OLD.py // \n')
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
            ###start by top nodes ...
            for j in range(0,Npttop):
                np=np+1
                geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_top['x'][j],DEM_top['y'][j])+r'}'+'; \n')
            ###Store the index of Domain top and bottom right corners
            DomainTopRightCornerPoint_Index = np
            DomainBottomRightCornerPoint_Index = np + 1
            ###...continu with bottom nodes
            for j in range(0,Nptbed):
                np=np+1
                geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_bed['x'][j],DEM_bed['y'][j])+r'}'+'; \n')
            ###Store the index of Domain bottom left corner
            DomainBottomLeftCornerPoint_Index = np
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ###...now write the coordinates (x,y,0,lc) FOR THE TUNNEL CONTOUR
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ####For the rectangle and half circle case, only a few points are required (it is important to turn counter-clokwise for the half-circle case)
            if Shape == 'Rectangle' or Shape == 'Circle':
                ##Write tunnel bottom left corner and store index
                np = np + 1
                TunnelBottomLeftCornerPoint_Index = np
                if AfterIncision:
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, BottomTunnel_altitude) + r'}' + '; \n')
                else:
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude) + r'}' + '; \n')
                ###Write tunnel bottom right corner and store index
                np = np + 1
                TunnelBottomRightCornerPoint_Index = np
                if AfterIncision:
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, BottomTunnel_altitude) + r'}' + '; \n')
                else:
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude) + r'}' + '; \n')
                ###Middle points on walls + Top corners
                if Shape == 'Rectangle':
                    ###Write middle point on right wall of tunnel and store index
                    np = np + 1
                    TunnelMiddlePointRightWall_Index = np
                    if AfterIncision:
                        if Incurved_Walls:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2 + xwallshift, (Floor_altitude+2.5-BottomTunnel_altitude)/2+BottomTunnel_altitude) + r'}' + '; \n')
                        else:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, (Floor_altitude+2.5-BottomTunnel_altitude)/2+BottomTunnel_altitude) + r'}' + '; \n')
                    else:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5/2) + r'}' + '; \n')
                    ###Write top right corner of tunnel and store index
                    np = np + 1
                    TunnelTopRightCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5) + r'}' + '; \n')
                    ###Write top left corner of tunnel and store index
                    np = np + 1
                    TunnelTopLeftCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5) + r'}' + '; \n')
                    ###Write middle point on left wall of tunnel and store index
                    np = np + 1
                    TunnelMiddlePointLeftWall_Index = np
                    if AfterIncision:
                        if Incurved_Walls:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2 - xwallshift, (Floor_altitude+2.5-BottomTunnel_altitude)/2+BottomTunnel_altitude) + r'}' + '; \n')
                        else:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, (Floor_altitude+2.5-BottomTunnel_altitude)/2+BottomTunnel_altitude) + r'}' + '; \n')
                    else:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5/2) + r'}' + '; \n')
                elif Shape == 'Circle':
                    ###Write middle point on right wall of tunnel and store index
                    np = np + 1
                    TunnelMiddlePointRightWall_Index = np
                    if AfterIncision:
                        if Incurved_Walls:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25 + Width / 2 + xwallshift,(Floor_altitude+2.5-0.5*Width - BottomTunnel_altitude) / 2 + BottomTunnel_altitude) + r'}' + '; \n')
                        else:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25 + Width / 2, (Floor_altitude+2.5-0.5*Width - BottomTunnel_altitude) / 2 + BottomTunnel_altitude) + r'}' + '; \n')
                    else:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25 + Width / 2, (Floor_altitude+2.5-0.5*Width - Floor_altitude) / 2 + Floor_altitude) + r'}' + '; \n')
                    ###Write top right corner of tunnel and store index
                    np = np + 1
                    TunnelTopRightCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
                    ###Write top left corner of tunnel and store index
                    np = np + 1
                    TunnelTopLeftCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
                    ###If circle, requires circle center as well
                    np = np + 1
                    TunnelRoofCenterPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25, Floor_altitude+2.5-0.5*Width) + r'}' + '; \n')
                    ###Write middle point on left wall of tunnel and store index
                    np = np + 1
                    TunnelMiddlePointLeftWall_Index = np
                    if AfterIncision:
                        if Incurved_Walls:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2 - xwallshift, (Floor_altitude+2.5-0.5*Width-BottomTunnel_altitude)/2 + BottomTunnel_altitude) + r'}' + '; \n')
                        else:
                            geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, (Floor_altitude+2.5-0.5*Width-BottomTunnel_altitude)/2 + BottomTunnel_altitude) + r'}' + '; \n')
                    else:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, (Floor_altitude+2.5-0.5*Width-Floor_altitude)/2 + Floor_altitude) + r'}' + '; \n')
            ####For the ovoid shape more points are required
            elif Shape == 'Ovoide':
                ####Get points coordinates from the DEM dataframe
                if AfterIncision: ###If tunnel after incision we must clean the bottom part of the ovoide which will have melted away after water incision
                    ###Find the altitude of the right/left extremities of the ovoide
                    ymin_ovoid=DEM_Ovoid[DEM_Ovoid['xt']==max(DEM_Ovoid['xt'])]['yt']
                    ###Keep only nodes of the ovoid that are above this altitude
                    DEM_Ovoid_Roof=DEM_Ovoid[DEM_Ovoid['yt'].values>=ymin_ovoid.values]
                    ###Reverse DEM_Ovoid_Roof dataframe to have decreasing x (counter-clockwise)
                    DEM_Ovoid_Roof = DEM_Ovoid_Roof[::-1]
                    DEM_Ovoid_Roof= DEM_Ovoid_Roof.reset_index(drop=True)
                    ##Write tunnel bottom left corner and store index
                    np = np + 1
                    TunnelBottomLeftCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(min(DEM_Ovoid['xt']), BottomTunnel_altitude) + r'}' + '; \n')
                    ###Write tunnel bottom right corner and store index
                    np = np + 1
                    TunnelBottomRightCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(max(DEM_Ovoid['xt']), BottomTunnel_altitude) + r'}' + '; \n')
                    ###Write middle point on right wall of tunnel and store index
                    np = np + 1
                    TunnelMiddlePointRightWall_Index = np
                    if Incurved_Walls:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(max(DEM_Ovoid['xt']) + xwallshift,(ymin_ovoid.values[0] - BottomTunnel_altitude) / 2 + BottomTunnel_altitude) + r'}' + '; \n')
                    else:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(max(DEM_Ovoid['xt']), (ymin_ovoid.values[0] - BottomTunnel_altitude) / 2 + BottomTunnel_altitude) + r'}' + '; \n')
                    ####Now write the coordinates of the roof point from the DEM dataframe
                    TunnelTopRightCornerPoint_Index = np + 1
                    for j in range(len(DEM_Ovoid_Roof)):
                        np = np + 1
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(DEM_Ovoid_Roof['xt'][j],DEM_Ovoid_Roof['yt'][j]) + r'}' + '; \n')
                    TunnelTopLeftCornerPoint_Index = np
                    ###Write middle point on left wall of tunnel and store index
                    np = np + 1
                    TunnelMiddlePointLeftWall_Index = np
                    if Incurved_Walls:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(min(DEM_Ovoid['xt']) - xwallshift,(ymin_ovoid.values[0] - BottomTunnel_altitude) / 2 + BottomTunnel_altitude) + r'}' + '; \n')
                    else:
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(min(DEM_Ovoid['xt']), (ymin_ovoid.values[0] - BottomTunnel_altitude) / 2 + BottomTunnel_altitude) + r'}' + '; \n')
                else: ###Here we consider the ovoid in its initial shape
                      ###It follows that middle points, bottom right and bottom left are not defined (bottom right/left actually corresponds to top right/left)
                    TunnelTopRightCornerPoint_Index = np+1
                    for j in range(len(DEM_Ovoid)):
                        np = np + 1
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(DEM_Ovoid['xt'][j], DEM_Ovoid['yt'][j]) + r'}' + '; \n')
                    TunnelTopLeftCornerPoint_Index = np-1 ##Because the last point of the intial ovoid is the floor middle point
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
                for j in range(DomainBottomRightCornerPoint_Index,DomainBottomLeftCornerPoint_Index):
                    geo.write('{},'.format(j))
                geo.write('{}'.format(DomainBottomLeftCornerPoint_Index) + '}; \n')
                ###The spline corresponding to right side of the domain
                geo.write('Spline(4) = {')
                geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
                geo.write('1}; \n')
                ###Now the spline for the tunnel (important to turn counter-clockwise for the half circle case)
                if Shape == 'Rectangle' or Shape == 'Circle':
                    ###Bottom line
                    geo.write('Spline(5) = {')
                    geo.write('{},'.format(TunnelBottomLeftCornerPoint_Index))
                    geo.write('{}'.format(TunnelBottomRightCornerPoint_Index))
                    geo.write('}; \n')
                    ###Right wall
                    geo.write('Spline(6) = {')
                    geo.write('{},'.format(TunnelBottomRightCornerPoint_Index))
                    geo.write('{},'.format(TunnelMiddlePointRightWall_Index))
                    geo.write('{}'.format(TunnelTopRightCornerPoint_Index))
                    geo.write('}; \n')
                    ###Top line (or half circle)
                    if Shape == 'Rectangle':
                        geo.write('Spline(7) = {')
                        geo.write('{},'.format(TunnelTopRightCornerPoint_Index))
                        geo.write('{}'.format(TunnelTopLeftCornerPoint_Index))
                        geo.write('}; \n')
                    elif Shape == 'Circle':
                        geo.write('Circle(7) = {')
                        geo.write('{},'.format(TunnelTopRightCornerPoint_Index))
                        geo.write('{},'.format(TunnelRoofCenterPoint_Index))
                        geo.write('{}'.format(TunnelTopLeftCornerPoint_Index))
                        geo.write('}; \n')
                    ###Left wall
                    geo.write('Spline(8) = {')
                    geo.write('{},'.format(TunnelTopLeftCornerPoint_Index))
                    geo.write('{},'.format(TunnelMiddlePointLeftWall_Index))
                    geo.write('{}'.format(TunnelBottomLeftCornerPoint_Index))
                    geo.write('}; \n')
                elif Shape == 'Ovoide':
                    if AfterIncision:
                        ###Bottom line
                        geo.write('Spline(5) = {')
                        geo.write('{},'.format(TunnelBottomLeftCornerPoint_Index))
                        geo.write('{}'.format(TunnelBottomRightCornerPoint_Index))
                        geo.write('}; \n')
                        ###Right wall
                        geo.write('Spline(6) = {')
                        geo.write('{},'.format(TunnelBottomRightCornerPoint_Index))
                        geo.write('{},'.format(TunnelMiddlePointRightWall_Index))
                        geo.write('{}'.format(TunnelTopRightCornerPoint_Index))
                        geo.write('}; \n')
                        ###Top line (or half circle)
                        geo.write('Spline(7) = {')
                        for j in range(TunnelTopRightCornerPoint_Index, TunnelTopLeftCornerPoint_Index):
                            geo.write('{},'.format(j))
                        geo.write('{}'.format(TunnelTopLeftCornerPoint_Index) + '}; \n')
                        ###Left wall
                        geo.write('Spline(8) = {')
                        geo.write('{},'.format(TunnelTopLeftCornerPoint_Index))
                        geo.write('{},'.format(TunnelMiddlePointLeftWall_Index))
                        geo.write('{}'.format(TunnelBottomLeftCornerPoint_Index))
                        geo.write('}; \n')
                    else:
                        ###Ovoid walls without floor
                        geo.write('Spline(5) = {')
                        for j in range(TunnelTopRightCornerPoint_Index, TunnelTopLeftCornerPoint_Index):
                            geo.write('{},'.format(j))
                        geo.write('{}'.format(TunnelTopLeftCornerPoint_Index) + '}; \n')
                        ###Ovoid floor
                        geo.write('Spline(6) = {')
                        geo.write('{},'.format(TunnelTopLeftCornerPoint_Index))
                        geo.write('{},'.format(TunnelTopLeftCornerPoint_Index+1)) ##The floor middle point
                        geo.write('{}'.format(TunnelTopRightCornerPoint_Index))
                        geo.write('}; \n')
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ###Now build the line loops
                geo.write('Line Loop(9) = {1, 2, 3, 4}; \n') ###domain contour
                if Shape == 'Rectangle' or Shape == 'Circle': ###tunnel contour
                    geo.write('Line Loop(10) = {5, 6, 7, 8}; \n')
                elif Shape == 'Ovoide':
                    if AfterIncision:
                        geo.write('Line Loop(10) = {5, 6, 7, 8}; \n')
                    else:
                        geo.write('Line Loop(10) = {5, 6}; \n')
                geo.write('Plane Surface(11) = {9,10}; \n')
                geo.write('Physical Line(12) = {1}; \n') ##BC top
                geo.write('Physical Line(13) = {2}; \n') ##BC right side
                geo.write('Physical Line(14) = {3}; \n') ##BC bottom
                geo.write('Physical Line(15) = {4}; \n') ##BC left side
                if Shape == 'Ovoide' and not(AfterIncision):
                    geo.write('Physical Line(16) = {5, 6}; \n') ##BC tunnel (Free surface
                else:
                    geo.write('Physical Line(16) = {5, 6, 7, 8}; \n') ##BC tunnel (Free surface)
                geo.write('Physical Surface(17) = {11}; \n')
            geo.close()

