# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 18:33:47 2016

@author: claire
"""

import os, sys, datetime, string
import numpy as np
#from netCDF4 import Dataset as NetCDFFile
import numpy.ma as ma
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from pylab import *
from scipy.io.netcdf import netcdf_file as Dataset 
import matplotlib as mpl

def makeMap(lonStart,lonEnd,latStart,latEnd,name,stLon,stLat,outfile):
    """Map the result on taiwan bagground map 
    *input :
    -lonStart, float ; long low limit of the map
    -lonEnd, float ; lat high limit of the map
    -latStart, float ; lat low limit of the map
    -latEnd, float ; lat high limit of the map
    -name; str,  name of the map
    -stLonList ; list of longitude of point to be represented
    -stLat, List ; list of latitude of point to be represented
    -outfile, str, outfile name and folder
    *output :
    draw map of the input point
    *exemple :
        lonStart= [119.5]
        lonEnd =[122.5]
        latStart =[21.5]
        latEnd = [26]
        name = ['Taiwan']
        Long56 = [119.7,120,122.5]
        Lat56 = [22.5,23,23.7]
        outfile = 'MAP_EQ_4_7distmax_200'
        makeMap(lonStart[0],lonEnd[0],latStart[0],latEnd[0],name[0],Long56,Lat56,outfile)
    """
    plt.figure(figsize=(10,10))
    """Get the etopo2 data"""
    etopo1name='/home/claire/PHD/Working/Data/Earthbathy/ETOPO1_Bed_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')
    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]
    
    res = findSubsetIndices(latStart-5,latEnd+5,lonStart-40,lonEnd+10,lats,lons)
    print 'longggggggg llaaaaat',lons[res[0]:res[1]],lats[res[2]:res[3]]
    lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    print "Extracted data for area %s : (%s,%s) to (%s,%s)"%(name,lon.min(),lat.min(),lon.max(),lat.max())
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)
    levels=[-6000,-5000,-3000, -2000, -1500, -1000,-500, -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, -25, -15, -10, -5,-2,-1, 0]
    if lonStart<0 and lonEnd<0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0
    map = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,llcrnrlon=lonStart,urcrnrlon=lonEnd,rsphere=(6378137.00,6356752.3142),resolution='h',area_thresh=1000.,projection='lcc',lat_1=latStart,lon_0=lon_0)
    x, y = map(lon,lat)
    map.drawcoastlines()
    map.drawcountries()
    map.fillcontinents(color='grey',alpha=0.3)
    map.drawmeridians(np.arange(lons.min(),lons.max(),2),labels=[0,0,0,1])
    map.drawparallels(np.arange(lats.min(),lats.max(),2),labels=[1,0,0,0])
    #map.bluemarble()
    
    CS1 = map.contourf(x,y,bathySmoothed,levels,
    cmap=LevelColormap(levels,cmap=cm.Blues_r),
    extend='upper',
    alpha=0.6,
    origin='lower')
    CS1.axis='tight'
    """Plot the station as a position dot on the map"""
    for i in xrange(len(stLat)):
        xpt,ypt = map(stLon[i],stLat[i])
        map.plot([xpt],[ypt],marker='*',color='r', markersize=10)
    #plt.text(xpt+100000,ypt+100000,name)
    plt.title('Area %s'%(name))
    plotfile='/home/claire/PHD/Working/%s.png'%outfile
    plt.savefig(plotfile,dpi=300,orientation='portrait')
    plt.show()
    return 
    
def makeMapValues(lonStart,lonEnd,latStart,latEnd,name,stLon,stLat,lvalues,lSize,legend,outfile):
    """Map the result on taiwan bagground map 
    *input :
    -lonStart, float ; long low limit of the map
    -lonEnd, float ; lat high limit of the map
    -latStart, float ; lat low limit of the map
    -latEnd, float ; lat high limit of the map
    -name; str,  name of the map
    -stLonList ; list of longitude of point to be represented
    -stLat, List ; list of latitude of point to be represented
    -lvalues, list, list of the values to be represented ex: Ml, Mean peak to Peak 
    -outfile, str, outfile name and folder
    *output :
    draw map of the input point
    *exemple :
        lonStart= [119.5]
        lonEnd =[122.5]
        latStart =[21.5]
        latEnd = [26]
        name = ['Taiwan']
        Long56 = [119.7,120,122.5]
        Lat56 = [22.5,23,23.7]
        outfile = 'MAP_EQ_4_7distmax_200'
        makeMap(lonStart[0],lonEnd[0],latStart[0],latEnd[0],name[0],Long56,Lat56,outfile)
    """
    fig=plt.figure(figsize=(10,10))
    """Get the etopo2 data"""
    etopo1name='/home/claire/PHD/Working/Data/Earthbathy/ETOPO1_Bed_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')
    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]
    
    res = findSubsetIndices(latStart-5,latEnd+5,lonStart-40,lonEnd+10,lats,lons)
    lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    print "Extracted data for area %s : (%s,%s) to (%s,%s)"%(name,lon.min(),lat.min(),lon.max(),lat.max())
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)
    levels=[-6000,-5000,-3000, -2000, -1500, -1000,-500, -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, -25, -15, -10, -5,-2,-1, 0]
    if lonStart<0 and lonEnd<0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0
    #Continent map and map border ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    map = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,llcrnrlon=lonStart,urcrnrlon=lonEnd,rsphere=(6378137.00,6356752.3142),resolution='h',area_thresh=1000.,projection='lcc',lat_1=latStart,lon_0=lon_0)
    x, y = map(lon,lat)
    map.drawcoastlines()
    map.drawcountries()
    map.fillcontinents(color='grey',alpha=0.3)
    map.drawmeridians(np.arange(lons.min(),lons.max(),2),labels=[0,0,0,1])
    map.drawparallels(np.arange(lats.min(),lats.max(),2),labels=[1,0,0,0])
    
    #Ocean map~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CS1 = map.contourf(x,y,bathySmoothed,levels,
    cmap=LevelColormap(levels,cmap=cm.Blues_r),
    extend='upper',
    alpha=0.6,
    origin='lower')
    #points to be represented ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CS1.axis='tight'
    """Plot the station as a position dot on the map"""
    Norm = MidpointNormalize(midpoint=1)
    xpt,ypt = map(stLon,stLat)
    cmap = plt.cm.get_cmap('seismic')
    
    pcm = map.scatter([xpt],[ypt],c=lvalues,norm=Norm,marker='o',s=lSize,edgecolors='none',cmap =cmap)
    xst,yst = map(121.38,23.69)
    map.plot(xst,yst,color ='k',marker='^',markersize=10)
    
    plt.title('%s'%(name))
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(pcm, cbar_ax)
    plt.title(legend,'right')
    #plt.text(xpt+100000,ypt+100000,name)
    plotfile='/home/claire/PHD/Working/%s.png'%outfile
    plt.savefig(plotfile,dpi=300,orientation='portrait')
    plt.show()
    return 
    
def makeMapMl(lonStart,lonEnd,latStart,latEnd,name,stLon,stLat,lMl,outfile):
    """Map the result magnitude of the EQ on taiwan bagground map 
    *input :
    -lonStart, float ; long low limit of the map
    -lonEnd, float ; lat high limit of the map
    -latStart, float ; lat low limit of the map
    -latEnd, float ; lat high limit of the map
    -name; str,  name of the map
    -stLonList ; list of longitude of point to be represented
    -stLat, List ; list of latitude of point to be represented
    -lvalues, list, list of the values to be represented ex: Ml, Mean peak to Peak 
    -outfile, str, outfile name and folder
    *output :
    draw map of the input point
    *exemple :
        lonStart= [119.5]
        lonEnd =[122.5]
        latStart =[21.5]
        latEnd = [26]
        name = ['Taiwan']
        Long56 = [119.7,120,122.5]
        Lat56 = [22.5,23,23.7]
        outfile = 'MAP_EQ_4_7distmax_200'
        makeMap(lonStart[0],lonEnd[0],latStart[0],latEnd[0],name[0],Long56,Lat56,outfile)
    """
    fig=plt.figure(figsize=(10,10))
    """Get the etopo2 data"""
    etopo1name='/home/claire/PHD/Working/Data/Earthbathy/ETOPO1_Bed_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')
    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]
    
    res = findSubsetIndices(latStart-5,latEnd+5,lonStart-40,lonEnd+10,lats,lons)
    print 'longggggggg llaaaaat',lons[res[0]:res[1]],lats[res[2]:res[3]]
    lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    print "Extracted data for area %s : (%s,%s) to (%s,%s)"%(name,lon.min(),lat.min(),lon.max(),lat.max())
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplace_filter(bathy,M=None)
    levels=[-6000,-5000,-3000, -2000, -1500, -1000,-500, -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, -25, -15, -10, -5,-2,-1, 0]
    if lonStart<0 and lonEnd<0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0
    #Continent map and map border ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    map = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,llcrnrlon=lonStart,urcrnrlon=lonEnd,rsphere=(6378137.00,6356752.3142),resolution='h',area_thresh=1000.,projection='lcc',lat_1=latStart,lon_0=lon_0)
    x, y = map(lon,lat)
    map.drawcoastlines()
    map.drawcountries()
    map.fillcontinents(color='grey',alpha=0.3)
    map.drawmeridians(np.arange(lons.min(),lons.max(),2),labels=[0,0,0,1])
    map.drawparallels(np.arange(lats.min(),lats.max(),2),labels=[1,0,0,0])
    
    #Ocean map~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CS1 = map.contourf(x,y,bathySmoothed,levels,
    cmap=LevelColormap(levels,cmap=cm.Blues_r),
    extend='upper',
    alpha=0.6,
    origin='lower')
    #points to be represented ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CS1.axis='tight'
    """Plot the station as a position dot on the map"""
    xpt,ypt = map(stLon,stLat)
    pcm = map.scatter([xpt],[ypt],color=lMl,marker='*',s=100,edgecolors='none')
    xst,yst = map(121.38,23.69)
    map.plot(xst,yst,color ='k',marker='^',markersize=10)
    plt.title('%s'%(name))
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cmap = mpl.colors.ListedColormap(['turquoise','mediumblue','limegreen','forestgreen' ,'gold' , 'darkgoldenrod','red' ,'darkred'])
    cmap.set_over('0.25')
    cmap.set_under('0.75')
    
    # If a ListedColormap is used, the length of the bounds array must be
    # one greater than the length of the color list.  The bounds must be
    # monotonically increasing.
    bounds = [3,3.5,4,4.5,5,5.5,6,6.5]
    norm = mpl.colors.Normalize(3,7)
    cb2 = mpl.colorbar.ColorbarBase(ax = cbar_ax,cmap=cmap,norm =norm,ticks=bounds)
    plt.title('Ml')
    #plt.text(xpt+100000,ypt+100000,name)
    
    plotfile='/home/claire/PHD/Working/%s.png'%outfile
    plt.savefig(plotfile,dpi=300,orientation='portrait')
    plt.show()
    return 
    
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
    

def findSubsetIndices(min_lat,max_lat,min_lon,max_lon,lats,lons):
    """Array to store the results returned from the function"""
    res=np.zeros((4),dtype=np.float64)
    minLon=min_lon; maxLon=max_lon
    distances1 = []; distances2 = []
    indices=[]; index=1
    for point in lats:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1
    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])
    distances1 = []; distances2 = []; index=1
        
    for point in lons:
        s1 = maxLon-point # (vector subtract)
        s2 = minLon-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1
    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])
    
    """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
    minJ=indices[1][2]
    maxJ=indices[0][2]
    minI=indices[3][2]
    maxI=indices[2][2]
    res[0]=minI; res[1]=maxI; res[2]=minJ; res[3]=maxJ;
    return res


def laplace_X(F,M):
    """1D Laplace Filter in X-direction (axis=1)"""
    jmax, imax = F.shape
    # Add strips of land
    F2 = np.zeros((jmax, imax+2), dtype=F.dtype)
    F2[:, 1:-1] = F
    M2 = np.zeros((jmax, imax+2), dtype=M.dtype)
    M2[:, 1:-1] = M
    MS = M2[:, 2:] + M2[:, :-2]
    FS = F2[:, 2:]*M2[:, 2:] + F2[:, :-2]*M2[:, :-2]
    return np.where(M > 0.5, (1-0.25*MS)*F + 0.25*FS, F)

def laplace_Y(F,M):
    """1D Laplace Filter in Y-direction (axis=1)"""
    jmax, imax = F.shape
    # Add strips of land
    F2 = np.zeros((jmax+2, imax), dtype=F.dtype)
    F2[1:-1, :] = F
    M2 = np.zeros((jmax+2, imax), dtype=M.dtype)
    M2[1:-1, :] = M
    MS = M2[2:, :] + M2[:-2, :]
    FS = F2[2:, :]*M2[2:, :] + F2[:-2, :]*M2[:-2, :]
    return np.where(M > 0.5, (1-0.25*MS)*F + 0.25*FS, F)

def laplace_filter(F, M=None):
    if M == None:
        M = np.ones_like(F)

    return 0.5*(laplace_X(laplace_Y(F, M), M)+laplace_Y(laplace_X(F, M), M))

# ------------------
# Plot land mask
# ------------------

def landmask(M, color='0.8'):

   # Make a constant colormap, default = grey
   constmap = pl.matplotlib.colors.ListedColormap([color])

   jmax, imax = M.shape
   # X and Y give the grid cell boundaries,
   # one more than number of grid cells + 1
   # half integers (grid cell centers are integers)
   X = -0.5 + pl.arange(imax+1)
   Y = -0.5 + pl.arange(jmax+1)

   # Draw the mask by pcolor
   M = ma.masked_where(M > 0, M)
   pl.pcolor(X, Y, M, shading='flat', cmap=constmap)

# -------------
# Colormap
# -------------

def LevelColormap(levels, cmap=None):
    """Make a colormap based on an increasing sequence of levels"""
   
    # Start with an existing colormap
    if cmap == None:
        cmap = pl.get_cmap()

    # Spread the colours maximally
    nlev = len(levels)
    S = pl.arange(nlev, dtype='float')/(nlev-1)
    A = cmap(S)

    # Normalize the levels to interval [0,1]
    levels = pl.array(levels, dtype='float')
    L = (levels-levels[0])/(levels[-1]-levels[0])

    # Make the colour dictionary
    R = [(L[i], A[i,0], A[i,0]) for i in xrange(nlev)]
    G = [(L[i], A[i,1], A[i,1]) for i in xrange(nlev)]
    B = [(L[i], A[i,2], A[i,2]) for i in xrange(nlev)]
    cdict = dict(red=tuple(R),green=tuple(G),blue=tuple(B))

    # Use
    return matplotlib.colors.LinearSegmentedColormap('%s_levels' % cmap.name, cdict, 256)
        

