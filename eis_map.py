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
outfilebase = 'eis_map_'
lts_outfilebase = 'lts_map_'
casename = sc

qs850 = filebase + casename+'/'+outfilebase+'_qs850.nc'
temp700 = filebase + casename + '/'+outfilebase+'_temp700.nc'
tempsurf = filebase + casename+'/'+outfilebase+'_tempsurf.nc'
tempsum = filebase + casename+'/'+outfilebase+'_tempsum.nc'
z700 = filebase + casename + '/'+outfilebase+'_z700.nc'

lcl = filebase + casename+'/'+outfilebase+'_lcl.nc'

# check directly if the mergetime file exists
if os.path.isdir(filebase+casename):
	if not os.path.isfile(lcl):
		infile = filebase + casename+'/merged.nc'

		if not os.path.isfile(qs850):
			# calc saturation mixing ratio at 850 hPa
			syscall = 'cdo timmean -selyear,21/40 -sellevidx,23 -select,name=smr -expr,\'smr=(Q/RELHUM)\'  '+infile+' '+qs850
			os.system(syscall)
		
		if not os.path.isfile(temp700):
			# get temp at 700hPa
			syscall = 'cdo timmean -selyear,21/40 -sellevidx,21 -select,name=T '+infile+' '+temp700
			os.system(syscall)

		if not os.path.isfile(tempsurf):
			# get temp at surface
			syscall = 'cdo timmean -selyear,21/40 -sellevidx,26 -select,name=T '+infile+' '+tempsurf
			os.system(syscall)

		if not os.path.isfile(z700):
			# get height of 700hPa
			syscall = 'cdo timmean -selyear,21/40 -sellevidx,21 -select,name=Z3 '+infile+' '+z700
			os.system(syscall)

		if not os.path.isfile(lcl):
			# calc lifting condensation level
			syscall = 'cdo  timmean -selyear,21/40 -sellevidx,26 -select,name=lcl -expr,\'lcl=((0.02+(TS-273.15)/5)*0.1*(1-RELHUM))\'  '+infile+' '+lcl
			os.system(syscall)

		if not os.path.isfile(tempsum):
			# temp sum 
			syscall = 'cdo add '+temp700+' '+tempsurf + ' '+tempsum
			os.system(syscall)

lcl = []
z700 = []
lts = []
qs850 = []
tempsum = []
pal = []
eis = []
lons = []
lats = []

	
# GET LCL
dsloc = filebase + casename+'/'+outfilebase+'_lcl.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	lcl = ds.variables['lcl'][:]
	ds.close() #close the file

#GET z700
dsloc = filebase + casename+'/'+outfilebase+'_z700.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	z700 = ds.variables['Z3'][:]
	ds.close() #close the file

#GET lts
dsloc = filebase + casename+'/'+lts_outfilebase+'sub.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	lts = ds.variables['lts'][:]
	ds.close() #close the file


#GET qs850
dsloc = filebase + casename+'/'+outfilebase+'_qs850.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	qs850 = ds.variables['smr'][:]
	ds.close() #close the file	

#GET tempsum
dsloc = filebase + casename+'/'+outfilebase+'_tempsum.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	tempsum = ds.variables['T'][:]
	lons = ds.variables['lon'][:]
	lats = ds.variables['lat'][:]
	ds.close() #close the file

	
pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	

eis = lts - pal*((z700/1000)-lcl)


#create plot
fig = plt.figure()

# setup the map
m = Basemap(lat_0=0,lon_0=0)
m.drawcoastlines()
m.drawcountries()
	
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)

# Plot map
cs = m.pcolormesh(lons,lats, np.squeeze(eis), cmap='PRGn_r', latlon='True', vmin=-40,vmax=40)

# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label('Potential Temperature (K)')	

plt.title('Global Average Estimated Inversion Strength: Solar Constant = '+sc)

plt.show()

fig.savefig("eis_map_"+sc+".pdf", bbox_inches='tight')

