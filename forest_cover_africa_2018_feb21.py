#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 11:15:49 2018

# Congo basin tree fraction using 2018 dataset with gain

@author: earjba
"""
import numpy as np
import importlib
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from netCDF4 import Dataset, num2date
from mpl_toolkits import basemap
from jpros import readfiles
from jpros import harmonised
importlib.reload(readfiles)
importlib.reload(harmonised)
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import numpy as np
import pandas as pd
import tifffile as tiff
import h5py
import glob
import math
import gdal
import iris
from netCDF4 import Dataset as NetCDFFile
from PIL import Image
from pyhdf.SD import SD, SDC
from mpl_toolkits import basemap



def get_coords(gt, width, height):
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*gt[5]
    resx = gt[1]
    maxx = gt[0] + width*gt[1] + height*gt[2]
    maxy = gt[3]
    resy = gt[5]
    lon = np.arange(minx, maxx, resx)
    lat = np.arange(miny, maxy, -resy)
    return(lat, lon)


def regrid_data(var, lat, lon, target_res=2):
    if lat[0] > lat[-1]:
        #flip lat and associated index
        lat = lat[::-1]
        var = var[:, :, ::-1, :]
    new_lat = np.arange(lat[0], lat[-1]+(abs(lat[1]-lat[0])), target_res)
    new_lon = np.arange(lon[0], lon[-1]+(abs(lat[1]-lat[0])), target_res)
    lon_sub, lat_sub = np.meshgrid(new_lon, new_lat)
    var_rescale = np.empty((var.shape[0], var.shape[1],
                            len(new_lat), len(new_lon)))
    for yr in range(var.shape[0]):
        for mn in range(12):
            var_rescale[yr, mn, :, :] = basemap.interp(var[yr, mn, :, :],
                                                       lon, lat,
                                                       lon_sub, lat_sub,
                                                       order=1)
    return(var_rescale, new_lat, new_lon)


def get_lat_bounds(array1d):
    div = abs(array1d[1] - array1d[0])
    if array1d[0] < 0:
        extra_val = array1d[0] - div
        bounds1d = np.concatenate(([extra_val], array1d))
    else:
        extra_val = array1d[-1] - div
        bounds1d = np.concatenate((array1d, [extra_val]))
    bounds2d = np.hstack((bounds1d[:-1, np.newaxis], bounds1d[1:, np.newaxis]))
    bounds2d = bounds2d.astype('float')
    return(bounds2d)


def get_lon_bounds(array1d):
    div = abs(array1d[1] - array1d[0])
    extra_val = array1d[-1] + div
    bounds1d = np.concatenate((array1d, [extra_val]))
    bounds2d = np.hstack((bounds1d[:-1, np.newaxis], bounds1d[1:, np.newaxis]))
    bounds2d = bounds2d.astype('float')
    return(bounds2d)


def minus180_to_plus180(var, lon):
    if len(var.shape) < 4:
        raise TypeError('Variable not in correct format')
    else:
        l = int(var.shape[-1]/2)
        temp1 = var[:, :, :, 0:l]
        temp2 = var[:, :, :, l:]
        new_var = np.concatenate((temp2, temp1), axis=3)
        new_lon = np.arange(-180, 180, (abs(lon[1]-lon[0])))
    return(new_var, new_lon)


path = ('/nfs/a68/gyjcab/datasets/lapse_data_harmonised/Jan_2018/mon_1.0deg/')
# read in surface air temperture 2001- 2018
one_deg_cube = path+'tas_airs_mon_1.0deg_2003_2018.nc'
path = ('/nfs/a68/gyjcab/datasets/lapse_data_harmonised/Jan_2018/mon_0.25deg/')
# read in surface albedo
pt25_cube = path+'sal_clara_mon_0.25deg_1982_2015_direct_from_netcdf.nc'
pt5_cube = ('/nfs/a68/gyjcab/datasets/lapse_data_harmonised/Jan_2018/Final/'
# read in Global Precipitation Climatology Centre monthly total of Precipitation
            '0.5deg/pr_gpcc_mon_0.5deg_1983_2018.nc')
regrid_cube = one_deg_cube
#%%
def get_forest_cover_2018(res=0.05):
    # read canopy cover data
    forest_2000_path = '/nfs/a68/gyjcab/datasets/GFC_Hansen/v1.6/treecover2000/'

    # read forest loss year data
    year_loss_path = '/nfs/a68/gyjcab/datasets/GFC_Hansen/v1.6/lossYear/'


    # read each tile and down-scale data over Central Africa
    scale = int(res/abs(0.00025))
    ydim = (int(40/res))
    xdim1 = (int(10/res))
    xdim2 = (int(40/res))
    vdat = np.empty((ydim, xdim1))
    hdat = np.empty((ydim, xdim2))
    j = 0
    # 0-40E, 20-20S
    for nx in np.arange(0,40,10):
        i = 0
        for ny in np.arange(20,-20,-10):
            xstr = str(nx).zfill(3)
            ystr = str(ny).zfill(2)

            # First get tree cover 2000
            fname = 'Hansen_GFC-2018-v1.6_treecover2000_' +\
                     ystr + 'N_' + xstr + 'E.tif'
            if ny < 0:
                ystr = str(abs(ny)).zfill(2)
                fname = 'Hansen_GFC-2018-v1.6_treecover2000_' +\
                     ystr + 'S_' + xstr + 'E.tif'
            print(fname)
            # open tiff
            ds = gdal.Open(forest_2000_path+fname)
            band = ds.GetRasterBand(1)
            treeFrac_2000 = band.ReadAsArray()
            width = ds.RasterXSize
            height = ds.RasterYSize
            gt = ds.GetGeoTransform()
            lat, lon = get_coords(gt, width, height)

            # Then get loss data
            fname = 'Hansen_GFC-2018-v1.6_lossyear_' +\
                     ystr + 'N_' + xstr + 'E.tif'
            if ny < 0:
                ystr = str(abs(ny)).zfill(2)
                fname = 'Hansen_GFC-2018-v1.6_lossyear_' +\
                     ystr + 'S_' + xstr + 'E.tif'
            print(fname)
            # open tiff
            ds = gdal.Open(year_loss_path+fname)
            band = ds.GetRasterBand(1)
            loss_data = band.ReadAsArray()


            # Get mask of forest lost by 2018
            # where loss_data has values greater than 16, set to 0
            loss_data[loss_data>19] = 0
            # for all years, set value to 1. loss is now just binary
            loss_data[(loss_data>=1)&(loss_data<=19)] = 1

            treeFrac_2018 = treeFrac_2000.copy()
            treeFrac_2018[loss_data==1] = 0


            # Find mean forest cover per 0.25 degree pixel
            xx, yy = lon[::scale],lat[::scale]
            fc16_data = np.zeros((xdim1, xdim1))
            for ix in range(len(xx)):
                temp_lon_ix = np.where((lon==xx[ix]))[0]
                for iy in range(len(yy)):
                    temp_lat_iy = np.where((lat==yy[iy]))[0]

                    # Mean tree cover in 2000
#                    temp = (treeFrac_2000[temp_lat_iy[0]:temp_lat_iy[0]+scale, :]
#                                         [:, temp_lon_ix[0]:temp_lon_ix[0]+scale])
#                    print(np.nanmean(temp))

                    # Mean tree cover in 2018
                    data_trim = (treeFrac_2018[temp_lat_iy[0]:temp_lat_iy[0]+scale, :]
                                              [:, temp_lon_ix[0]:temp_lon_ix[0]+scale])

                    mean_tree_cover = np.nanmean(data_trim)
#                    print(mean_tree_cover)
                    fc16_data[iy, ix] = mean_tree_cover
#                assert False
            vdat[i*xdim1:(i+1)*xdim1,:] = fc16_data
            i += 1
        hdat[:,j*xdim1:(j+1)*xdim1] = vdat
        j += 1
    lat = np.arange(20, -20, -res)
    lon = np.arange(0, 40, res)
    return(hdat, lat, lon)

# LANDSAT forest loss
res = 0.05
hdat, lat, lon = get_forest_cover_2018(res=res)
plt.imshow(hdat)

if res == 0.25:
    regrid_cube = pt25_cube
elif res == 0.05:
    regrid_cube = ('/nfs/see-fs-02_users/earjba/python_scripts/'
                   'deforestation_analysis/temp_cube_5km.nc')

args = [hdat, lat, lon, 2018, 0, '1']
kwargs = {'standard_name': 'area_fraction',
          'long_name': 'Tree Cover Fraction (year 2018)',
          'short_name': 'treeFrac2018_' + str(res),
          'product': 'hansen-landsat',
          'regrid': True,
          'latlon_bounds': True,
          'time_bounds': False,
          'regrid_cube': regrid_cube}
harmonised.write_netcdf(*args, **kwargs)


### fix artefacts
path = '/nfs/a68/ee13c2s/python/treefrac_africa/2018/feb21/'
if res == 0.05:
    fname = 'treeFrac2018_0.05_hansen-landsat_mon_0.05deg.nc'
if res == 0.25:
    fname = 'treeFrac2018_0.25_hansen-landsat_mon_0.25deg.nc'
cube = iris.load_cube(path + fname)
temp_lat = cube.coord('latitude').points
temp_lon = cube.coord('longitude').points
i, = np.where((temp_lat>lat.max())|(temp_lat<lat.min()))
j, = np.where((temp_lon>lon.max())|(temp_lon<lon.min()))
cube.data[i, :] = np.nan
cube.data[:, j] = np.nan
plt.imshow(cube.data, origin='lower')
plt.colorbar()
iris.save(cube, path+fname)
