"""
@author: jbrondex

Description:
------------
This file aims at creating a geo (gmsh input file) from two 2D (x,y) DEM files (top and bed) corresponding to 7 transects extracted from 3D DEMs of surface and bed of Grande Motte Glacier (Tignes)
using QGIS. For each transect we add a soupirail footprint considering 3 possible shapes for the initial tunnel: rectangle, rectangle + half-circle, ovoide. The soupirail is the shape of the tunnel after
incision of ice by water down to the bedrock when only a small space above the bedrock as survived the ice creep.

"""

# -*- coding: utf-8 -*-
# Create a geo (gmsh input file) file from two contour files
# the contour files contain the (x,y) coordinates of the ordered
# points defining the contour of the glacier and the contour of the cavity
#
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import interpolate
from scipy.interpolate import CubicSpline

###We define a function to interpolate bottom corners of soupirail on bed DEM
def Interpolate_on_bed(DEM_bed, xcorner):
    xbed = DEM_bed['x']
    xbed = xbed[::-1]  ###the scipy interpolation function requires the x to be sorted in increasing order
    ybed = DEM_bed['y']
    ybed = ybed[::-1]  ###so that xbed and ybed are consistent
    interpolation = interpolate.CubicSpline(xbed, ybed)
    ycorner = interpolation(xcorner)
    return ycorner

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
                file_name_output = 'GrandeMotte_Soupirail_{}_{}.geo'.format(Transect, Shape)
            else:
                file_name_output='GrandeMotte_Soupirail_{}_{}_W{}cm.geo'.format(Transect,Shape,Width_Name)
            geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
            geo.write('// This a a geo file created using the python script Makegeo_GrandeMotteSoupirail_Transect.py // \n')
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
            ###Store the index of Domain top and bottom right corner
            DomainTopRightCornerPoint_Index = np
            DomainBottomRightCornerPoint_Index = np + 1
            ###...continu with bottom nodes on the right side of soupirail
            ###First bottom right corner of tunnel must be projected on BedDEM
            if Shape == 'Rectangle' or Shape == 'Circle':
                xbottomrightcorner = 25+Width/2
            elif Shape == 'Ovoide':
                xbottomrightcorner = max(DEM_Ovoid['xt'])
            ybottomrightcorner = Interpolate_on_bed(DEM_bed,xbottomrightcorner)
            ###Same thing for bottom left corner
            if Shape == 'Rectangle' or Shape == 'Circle':
                xbottomleftcorner = 25 - Width / 2
            elif Shape == 'Ovoide':
                xbottomleftcorner = min(DEM_Ovoid['xt'])
            ybottomleftcorner = Interpolate_on_bed(DEM_bed, xbottomleftcorner)
            ###Then we must split the DEM_bed to keep only the part located on the right side of the soupirail
            DEM_bed_right = DEM_bed[DEM_bed['x']>xbottomrightcorner]
            Nptbed_right = len(DEM_bed_right.index)
            ###Write the points corresponding to the bed on the right of the soupirail
            for j in range(0,Nptbed_right):
                np=np+1
                geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lc'.format(DEM_bed_right['x'][j],DEM_bed_right['y'][j])+r'}'+'; \n')
            ###Now write the point corresponding to the bottom right corner of the soupirail and store index
            np = np + 1
            ChannelBottomRightCornerPoint_Index = np
            geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lcc'.format(xbottomrightcorner,ybottomrightcorner)+r'}'+'; \n')
            ###Now write the points required for the roof of the soupirail
            ####For the rectangle and circle case, only a few points are required
            if Shape == 'Rectangle': ###For rectangle we arbitrarily set height of walls equals to width
                ###Top right corner
                np = np + 1
                ChannelRoofRightCornerPoint_Index = np
                geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25+Width/2, ybottomrightcorner + Width) + r'}' + '; \n')
                ### Top Left corner
                np = np + 1
                ChannelRoofLeftCornerPoint_Index = np
                geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25-Width/2, ybottomleftcorner + Width) + r'}' + '; \n')
                ###Now write the point corresponding to the bottom left corner and store index
                np = np + 1
                ChannelBottomLeftCornerPoint_Index = np
                geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lcc'.format(xbottomleftcorner,ybottomleftcorner)+r'}'+'; \n')
            elif Shape == 'Circle': ###For circle there are no lateral walls
                ### Top right is also bottom right and has already been written
                ### Top Left is bottom left and needs to be written
                np = np + 1
                ChannelBottomLeftCornerPoint_Index = np
                geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(xbottomleftcorner,ybottomleftcorner) + r'}' + '; \n')
                ###For circle, requires circle center as well
                np = np + 1
                ChannelBottomCenterPoint_Index = np
                geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(25, 0.5*(ybottomrightcorner+ybottomleftcorner)) + r'}' + '; \n')
            ####For the ovoid shape more points are required
            elif Shape == 'Ovoide':
                ###Find the altitude of the right/left extremities of the ovoide
                ymin_ovoid=DEM_Ovoid[DEM_Ovoid['xt']==max(DEM_Ovoid['xt'])]['yt']
                ###Keep only nodes of the ovoid that are above this altitude
                DEM_Ovoid_Roof=DEM_Ovoid[DEM_Ovoid['yt'].values>=ymin_ovoid.values]
                ###Find which of the right or left corner of soupirail is highest (depends on local bed slope)
                IsRightHighest = ybottomrightcorner > ybottomleftcorner
                if IsRightHighest:
                    ### if right corner is highest then bring back all y of ovoid down to right corner ref
                    DEM_Ovoid_Roof_Corr = pd.DataFrame({'xr': DEM_Ovoid_Roof['xt'], 'yr': DEM_Ovoid_Roof['yt']-(ymin_ovoid.values-ybottomrightcorner)})
                    ### Top right is also bottom right and has already been written
                    ###Reverse DEM_Ovoid_Roof dataframe to have decreasing x
                    DEM_Ovoid_Roof_Corr = DEM_Ovoid_Roof_Corr[::-1]
                    DEM_Ovoid_Roof_Corr= DEM_Ovoid_Roof_Corr.reset_index(drop=True)
                    ####Now write the coordinates of the roof point from the DEM dataframe
                    for j in range(len(DEM_Ovoid_Roof_Corr)):
                        np = np + 1
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(DEM_Ovoid_Roof_Corr['xr'][j], DEM_Ovoid_Roof_Corr['yr'][j]) + r'}' + '; \n')
                    ChannelRoofLeftCornerPoint_Index = np
                    ###Now write the point corresponding to the bottom left corner and store index
                    np = np + 1
                    ChannelBottomLeftCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np)+r'{'+' {0}, {1}, 0.0, lcc'.format(xbottomleftcorner,ybottomleftcorner)+r'}'+'; \n')
                else: ##if left corner is the highest then bring back all y of ovoid down to left corner ref
                    DEM_Ovoid_Roof_Corr = pd.DataFrame({'xr': DEM_Ovoid_Roof['xt'], 'yr': DEM_Ovoid_Roof['yt'] - (ymin_ovoid.values - ybottomleftcorner)})
                    ### Top right is not bottom right and needs to be written
                    np = np + 1
                    ChannelRoofRightCornerPoint_Index = np
                    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(xbottomrightcorner, ybottomleftcorner) + r'}' + '; \n')
                    ###Reverse DEM_Ovoid_Roof dataframe to have decreasing x
                    DEM_Ovoid_Roof_Corr = DEM_Ovoid_Roof_Corr[::-1]
                    DEM_Ovoid_Roof_Corr= DEM_Ovoid_Roof_Corr.reset_index(drop=True)
                    ####Now write the coordinates of the roof point from the DEM dataframe
                    for j in range(len(DEM_Ovoid_Roof_Corr)):
                        np = np + 1
                        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(DEM_Ovoid_Roof_Corr['xr'][j], DEM_Ovoid_Roof_Corr['yr'][j]) + r'}' + '; \n')
                    ###In that case Top left is also bottom left and has already been written
                    ChannelBottomLeftCornerPoint_Index = np
            ###...continu with bottom nodes on the left side of channe
            ###Then we must split the DEM_bed to keep only the part located on the left side of the channel
            DEM_bed_left = DEM_bed[DEM_bed['x'] < xbottomleftcorner]
            Nptbed_left = len(DEM_bed_left.index)
            DEM_bed_left= DEM_bed_left.reset_index(drop=True)
            ###Write the points corresponding to the bed on the right of the channel
            for j in range(0, Nptbed_left):
                np = np + 1
                geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(DEM_bed_left['x'][j],DEM_bed_left['y'][j]) + r'}' + '; \n')
            DomainBottomLeftCornerPoint_Index = np

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
                ###The spline corresponding to the bed on the right side of the channel
                geo.write('Spline(3) = {')
                for j in range(DomainBottomRightCornerPoint_Index,ChannelBottomRightCornerPoint_Index):
                    geo.write('{},'.format(j))
                geo.write('{}'.format(ChannelBottomRightCornerPoint_Index) + '}; \n')
                ###The following depends on considered shape
                if Shape == 'Rectangle': ###if shape = rectangles there are walls on both sides
                    ###The spline corresponding to the channel right wall
                    geo.write('Spline(4) = {')
                    geo.write('{},'.format(ChannelBottomRightCornerPoint_Index))
                    geo.write('{}'.format(ChannelRoofRightCornerPoint_Index)+'}; \n')
                    ###The spline corresponding to the channel roof
                    geo.write('Spline(5) = {')
                    geo.write('{},'.format(ChannelRoofRightCornerPoint_Index))
                    geo.write('{}'.format(ChannelRoofLeftCornerPoint_Index))
                    geo.write('}; \n')
                    ###The spline corresponding to the channel left wall
                    geo.write('Spline(6) = {')
                    geo.write('{},'.format(ChannelRoofLeftCornerPoint_Index))
                    geo.write('{}'.format(ChannelBottomLeftCornerPoint_Index)+'}; \n')
                    ###The spline corresponding to the bed on the left side of the channel
                    geo.write('Spline(7) = {')
                    for j in range(ChannelBottomLeftCornerPoint_Index, DomainBottomLeftCornerPoint_Index):
                        geo.write('{},'.format(j))
                    geo.write('{}'.format(DomainBottomLeftCornerPoint_Index) + '}; \n')
                    ###The spline corresponding to left side of the domain
                    geo.write('Spline(8) = {')
                    geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
                    geo.write('1}; \n')
                    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    ###Now build the line loops
                    geo.write('Line Loop(9) = {1, 2, 3, 4, 5, 6, 7, 8}; \n')  ###domain contour
                    geo.write('Plane Surface(10) = {9}; \n')
                    geo.write('Physical Line(11) = {1}; \n')  ##BC top
                    geo.write('Physical Line(12) = {2}; \n')  ##BC domain right side
                    geo.write('Physical Line(13) = {3, 7}; \n')  ##BC bottom
                    geo.write('Physical Line(14) = {8}; \n')  ##BC domain left side
                    geo.write('Physical Line(15) = {4, 5, 6}; \n')  ##BC channel
                    geo.write('Physical Surface(16) = {10}; \n')
                elif Shape == 'Circle': ###if shape = circle there are no walls and top right/left are actually bottom right/left
                    geo.write('Circle(4) = {')
                    geo.write('{},'.format(ChannelBottomRightCornerPoint_Index))
                    geo.write('{},'.format(ChannelBottomCenterPoint_Index))
                    geo.write('{}'.format(ChannelBottomLeftCornerPoint_Index))
                    geo.write('}; \n')
                    ###The spline corresponding to the bed on the left side of the channel
                    geo.write('Spline(5) = {')
                    geo.write('{},'.format(ChannelBottomLeftCornerPoint_Index))
                    ### be carefull that point with index ChannelBottomLeftCornerPoint_Index+1 is center of circle so skip it
                    for j in range(ChannelBottomLeftCornerPoint_Index+2, DomainBottomLeftCornerPoint_Index):
                        geo.write('{},'.format(j))
                    geo.write('{}'.format(DomainBottomLeftCornerPoint_Index) + '}; \n')
                    ###The spline corresponding to left side of the domain
                    geo.write('Spline(6) = {')
                    geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
                    geo.write('1}; \n')
                    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    ###Now build the line loops
                    geo.write('Line Loop(7) = {1, 2, 3, 4, 5, 6}; \n')  ###domain contour
                    geo.write('Plane Surface(8) = {7}; \n')
                    geo.write('Physical Line(9) = {1}; \n')  ##BC top
                    geo.write('Physical Line(10) = {2}; \n')  ##BC domain right side
                    geo.write('Physical Line(11) = {3, 5}; \n')  ##BC bottom
                    geo.write('Physical Line(12) = {6}; \n')  ##BC domain left side
                    geo.write('Physical Line(13) = {4}; \n')  ##BC channel
                    geo.write('Physical Surface(14) = {8}; \n')
                elif Shape == 'Ovoide': ###if shape = ovoid, it depends of the slope of the bed
                    if IsRightHighest: ##no wall on the right but a wall on the left
                        geo.write('Spline(4) = {')
                        for j in range(ChannelBottomRightCornerPoint_Index,ChannelRoofLeftCornerPoint_Index):
                            geo.write('{},'.format(j))
                        geo.write('{}'.format(ChannelRoofLeftCornerPoint_Index) + '}; \n')
                        ###The spline corresponding to the channel left wall
                        geo.write('Spline(5) = {')
                        geo.write('{},'.format(ChannelRoofLeftCornerPoint_Index))
                        geo.write('{}'.format(ChannelBottomLeftCornerPoint_Index)+'}; \n')
                    else: ##a wall on the right but no wall on the left
                        ###The spline corresponding to the channel right wall
                        geo.write('Spline(4) = {')
                        geo.write('{},'.format(ChannelBottomRightCornerPoint_Index))
                        geo.write('{}'.format(ChannelRoofRightCornerPoint_Index)+'}; \n')
                        ###The spline corresponding to the roof
                        geo.write('Spline(5) = {')
                        for j in range(ChannelRoofRightCornerPoint_Index,ChannelBottomLeftCornerPoint_Index):
                            geo.write('{},'.format(j))
                        geo.write('{}'.format(ChannelBottomLeftCornerPoint_Index) + '}; \n')
                    ###The spline corresponding to the bed on the left side of the channel
                    geo.write('Spline(6) = {')
                    for j in range(ChannelBottomLeftCornerPoint_Index, DomainBottomLeftCornerPoint_Index):
                        geo.write('{},'.format(j))
                    geo.write('{}'.format(DomainBottomLeftCornerPoint_Index) + '}; \n')
                    ###The spline corresponding to left side of the domain
                    geo.write('Spline(7) = {')
                    geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
                    geo.write('1}; \n')
                    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    ###Now build the line loops
                    geo.write('Line Loop(8) = {1, 2, 3, 4, 5, 6, 7}; \n')  ###domain contour
                    geo.write('Plane Surface(9) = {8}; \n')
                    geo.write('Physical Line(10) = {1}; \n')  ##BC top
                    geo.write('Physical Line(11) = {2}; \n')  ##BC domain right side
                    geo.write('Physical Line(12) = {3, 6}; \n')  ##BC bottom
                    geo.write('Physical Line(13) = {7}; \n')  ##BC domain left side
                    geo.write('Physical Line(14) = {4, 5}; \n')  ##BC channel
                    geo.write('Physical Surface(15) = {9}; \n')
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            geo.close()

