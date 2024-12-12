"""
@author: jbrondex

Description:
------------
This file aims at creating a geo (gmsh input file) for a synthetic mesh corresponding to a symmetric half-domain. This will be
used to check how a non-penetration condition on tunnel/channel wall affects the closure of the tunnel/channel

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
#BottomTunnel_altitudes=[2792.2,2791.7,2791.2,2790.7,2790.2,2789.82] ### The elevation of tunnel bottom if Tunnel_After_Incision = True
BottomTunnel_altitudes=[2785.5,2785,2784.5,2784,2783.5,2783.12] ###case with incision to ~few cm above bed at transect 1 and then 1% slope
#### Possible shapes of the tunnel : Ovoid only for this case
Shape = 'Ovoide'
#### Do we want the half tunnel case or the half channel case ?
Case = 'Tunnel' ##'Tunnel' Or 'Channel'
###Build the ovoid DEM from parameterized curve method for the ovoid shape (based on RTM ppt)
t = np.linspace(-5, 5, 301)
x = 3.2 * t / (1 + t ** 2) ** 2
y = 3.2 / (1 + t ** 2) ** 2
df = pd.DataFrame({'xo': x, 'yo': y})
#### Remove points that are below 2.5m of maximum height
df = df[df['yo'] > np.max(df['yo']) - 2.5]
df = df[df['xo']>=0]
#### Keep only the right half of the shape
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
    DEM_bed = pd.read_csv(pathroot_mycode.joinpath(file_name_bed), names=Col_Names, delim_whitespace=True)
    ###Get averaged altitude of bed and of surface
    ymean_top=DEM_top['y'].mean()
    ymean_bed=DEM_bed['y'].mean()

    ######################################################
    ########## START WRITING THE .geo FILE  ##############
    ######################################################
    ###For the ovoid case, the tunnel must be centered based on floor altitude
    DEM_Ovoid = pd.DataFrame({'xt': df['xo'], 'yt': Floor_altitude - Floor_Ovoid + df['yo']})
    ####Get points coordinates from the DEM dataframe
    ###Find the altitude of the right/left extremities of the ovoide
    ymin_ovoid = DEM_Ovoid[DEM_Ovoid['xt'] == max(DEM_Ovoid['xt'])]['yt']
    ###Keep only nodes of the ovoid that are above this altitude
    DEM_Ovoid_Roof = DEM_Ovoid[DEM_Ovoid['yt'].values >= ymin_ovoid.values]
    ###Reverse DEM_Ovoid_Roof dataframe to have decreasing x (counter-clockwise)
    DEM_Ovoid_Roof = DEM_Ovoid_Roof[::-1]
    DEM_Ovoid_Roof = DEM_Ovoid_Roof.reset_index(drop=True)
    # Open the output file and write some general info
    file_name_output = 'GrandeMotte_Half{}SyntheSym_Deep_{}_{}.geo'.format(Case, Transect, Shape)
    geo = open(pathroot_mycode.joinpath(file_name_output), 'w')
    geo.write('// This a a geo file created using the python script Makegeo_GrandeMotte_SyntheticTunnelSym_Transect.py // \n')
    geo.write('Mesh.Algorithm=5; \n')
    geo.write('// To control the element size, one can directly modify the lc values in the geo file // \n')
    geo.write('lc = {} ; \n'.format(el_size))
    geo.write('lcc = {} ; \n'.format(el_sizec))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # write the points coordinates (x,y,0,lc) FOR THE DOMAIN CONTOUR
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###DEM top and bed have to be treated separately to have their own BCs for the domain contour
    np=1
    ###start by top nodes ...
    DomainTopLeftCornerPoint_Index = np
    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(0.0, ymean_top) + r'}' + '; \n')
    np = np + 1
    DomainTopRightCornerPoint_Index = np
    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(25.0, ymean_top) + r'}' + '; \n')
    ###... continue with bottom nodes
    np = np + 1
    DomainBottomRightCornerPoint_Index = np
    geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(25.0, ymean_bed) + r'}' + '; \n')
    np = np + 1
    DomainBottomLeftCornerPoint_Index = np
    if Case == 'Tunnel': ###If tunnel, bottom left of domain is on symmetry axis
        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lc'.format(0.0, ymean_bed) + r'}' + '; \n')
    elif Case == 'Channel': ###If Channel, bottom left of domain is also bottom left of tunnel
        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(max(DEM_Ovoid['xt']), ymean_bed) + r'}' + '; \n')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###...now write the coordinates (x,y,0,lc) FOR THE TUNNEL CONTOUR
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##Write tunnel bottom left corner and store index
    if Case == 'Tunnel': ###If tunnel, bottom left and bottom right based on floor altitude are required
        np = np + 1
        TunnelBottomLeftCornerPoint_Index = np
        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(0.0, BottomTunnel_altitude) + r'}' + '; \n')
        ###Write tunnel bottom right corner and store index
        np = np + 1
        TunnelBottomRightCornerPoint_Index = np
        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(max(DEM_Ovoid['xt']), BottomTunnel_altitude) + r'}' + '; \n')
    ####Now write the coordinates of the roof point from the DEM dataframe (same for tunnel and channel)
    TunnelTopRightCornerPoint_Index = np + 1
    ChannelTopRightCornerPoint_Index = np + 1
    for j in range(len(DEM_Ovoid_Roof)):
        np = np + 1
        geo.write('Point({}) = '.format(np) + r'{' + ' {0}, {1}, 0.0, lcc'.format(DEM_Ovoid_Roof['xt'][j],DEM_Ovoid_Roof['yt'][j]) + r'}' + '; \n')
    TunnelTopLeftCornerPoint_Index = np
    ChannelTopLeftCornerPoint_Index = np
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###Now build the lines from the points
    # if spline
    if spline:
        ###The spline corresponding to the top surface of the domain (same for channel and tunnel)
        geo.write('Spline(1) = {')
        geo.write('{},'.format(DomainTopLeftCornerPoint_Index))
        geo.write('{}'.format(DomainTopRightCornerPoint_Index)+'}; \n')
        ###The spline corresponding to right side of the domain (same for channel and tunnel)
        geo.write('Spline(2) = {')
        geo.write('{},'.format(DomainTopRightCornerPoint_Index))
        geo.write('{}'.format(DomainBottomRightCornerPoint_Index)+'}; \n')
        ###The spline corresponding to the bed of the domain (same for channel and tunnel although not same length)
        geo.write('Spline(3) = {')
        geo.write('{},'.format(DomainBottomRightCornerPoint_Index))
        geo.write('{}'.format(DomainBottomLeftCornerPoint_Index)+'}; \n')
        ###From there tunnel case requires more lines than channel case
        if Case == 'Tunnel':
            ###The spline corresponding to the bottom part of the symmetry axis
            geo.write('Spline(4) = {')
            geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
            geo.write('{}'.format(TunnelBottomLeftCornerPoint_Index)+'}; \n')
            ###The spline corresponding to the bottom part of the tunnel
            geo.write('Spline(5) = {')
            geo.write('{},'.format(TunnelBottomLeftCornerPoint_Index))
            geo.write('{}'.format(TunnelBottomRightCornerPoint_Index)+'}; \n')
            ###The spline corresponding to the right wall of the tunnel
            geo.write('Spline(6) = {')
            geo.write('{},'.format(TunnelBottomRightCornerPoint_Index))
            geo.write('{}'.format(TunnelTopRightCornerPoint_Index)+'}; \n')
            ###The spline corresponding to the roof of the tunnel
            geo.write('Spline(7) = {')
            for j in range(TunnelTopRightCornerPoint_Index, TunnelTopLeftCornerPoint_Index):
                geo.write('{},'.format(j))
            geo.write('{}'.format(TunnelTopLeftCornerPoint_Index) + '}; \n')
            ###The spline corresponding to the upper part of the symmetry axis
            geo.write('Spline(8) = {')
            geo.write('{},'.format(TunnelTopLeftCornerPoint_Index))
            geo.write('{}'.format(DomainTopLeftCornerPoint_Index)+'}; \n')
        elif Case == 'Channel':
            ###The spline corresponding to the right wall of the channel
            geo.write('Spline(4) = {')
            geo.write('{},'.format(DomainBottomLeftCornerPoint_Index))
            geo.write('{}'.format(ChannelTopRightCornerPoint_Index)+'}; \n')
            ###The spline corresponding to the roof of the tunnel
            geo.write('Spline(5) = {')
            for j in range(ChannelTopRightCornerPoint_Index, ChannelTopLeftCornerPoint_Index):
                geo.write('{},'.format(j))
            geo.write('{}'.format(ChannelTopLeftCornerPoint_Index) + '}; \n')
            ###The spline corresponding to the upper part of the symmetry axis
            geo.write('Spline(6) = {')
            geo.write('{},'.format(ChannelTopLeftCornerPoint_Index))
            geo.write('{}'.format(DomainTopLeftCornerPoint_Index)+'}; \n')
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ###Now build the line loops
        if Case == 'Tunnel':
            geo.write('Line Loop(9) = {1, 2, 3, 4, 5, 6, 7, 8}; \n') ###domain contour
            geo.write('Plane Surface(10) = {9}; \n')
            geo.write('Physical Line(11) = {1}; \n') ##BC top
            geo.write('Physical Line(12) = {2}; \n') ##BC right side
            geo.write('Physical Line(13) = {3}; \n') ##BC bottom
            geo.write('Physical Line(14) = {4, 8}; \n') ##BC symmetry axis
            geo.write('Physical Line(15) = {5, 6, 7}; \n') ##BC tunnel
            geo.write('Physical Surface(16) = {10}; \n')
        elif Case == 'Channel':
            geo.write('Line Loop(7) = {1, 2, 3, 4, 5, 6}; \n') ###domain contour
            geo.write('Plane Surface(8) = {7}; \n')
            geo.write('Physical Line(9) = {1}; \n') ##BC top
            geo.write('Physical Line(10) = {2}; \n') ##BC right side
            geo.write('Physical Line(11) = {3}; \n') ##BC bottom
            geo.write('Physical Line(12) = {6}; \n') ##BC symmetry axis
            geo.write('Physical Line(13) = {4,5}; \n') ##BC tunnel
            geo.write('Physical Surface(14) = {8}; \n')
    geo.close()

