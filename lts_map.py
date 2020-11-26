#!/usr/bin/env python3

import os
import sys
import numpy as np
import netCDF4
import operator
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

sc = sys.argv[1]

filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
outfilebase = 'lts_map_'
casename = sc

# calc temp
T700 = filebase + casename+'/'+outfilebase+'_700Temp.nc'
outfilesub = filebase + casename+'/'+outfilebase+'sub.nc'
outfile700 = filebase + casename+'/'+outfilebase+'700.nc'
outfile1000 = filebase + casename+'/'+outfilebase+'1000.nc'

# check directly if the lts file exists
if not os.path.isfile(outfilesub):
	infile = filebase + casename+'/merged.nc'

	# calc T at 700 hPa
	syscall = 'cdo select,name=T,PS -sellevel,696.79629 '+infile+ ' '+T700
	os.system(syscall)

	# calc potential temp 700
	syscall = 'cdo timmean -selyear,21/40 -select,name=lts -expr,\'lts=(T*(1000/696.79629)^0.286)\'  '+T700+' '+outfile700
	os.system(syscall)

	# calc potential temp at surface
	syscall = 'cdo timmean -selyear,21/40 -select,name=lts -expr,\'lts=(TS*(1000/(PS*0.01))^0.286)\'  '+infile+' '+outfile1000
	os.system(syscall)

	# difference = lower tropospheric stability
	syscall = 'cdo sub '+outfile700+' '+outfile1000 + ' '+outfilesub
	os.system(syscall)

#plot the data
if os.path.isfile(outfilesub):
	# open the file and get out the variable
	dsyear = netCDF4.Dataset(outfilesub)
	lons = dsyear.variables['lon'][:]
	lats = dsyear.variables['lat'][:]
	yr_lts = dsyear.variables['lts'][:]

	dsyear.close() #close the file

	#create plot
	fig = plt.figure()

	# setup the map
	m = Basemap(lat_0=0,lon_0=0)
	m.drawcoastlines()
	m.drawcountries()
	
	# Create 2D lat/lon arrays for Basemap
	lon2d, lat2d = np.meshgrid(lons, lats)

	# Plot of alb fraction with 11 contour intervals
	cs = m.pcolormesh(lons,lats,np.squeeze(yr_lts), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40)

	# This is the fix for the white lines between contour levels
	cs.set_edgecolor("face")

	# Add Colorbar
	cbar = m.colorbar(cs, location='bottom', pad="10%")
	cbar.set_label('Potential Temperature (K)')	

	plt.title('Global Average Lower Tropospheric Stability: '+ sc)
	plt.show()

	fig.savefig("lts_map_"+sc+".pdf", bbox_inches='tight')

