#!/usr/bin/env python3
"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""
import matplotlib as mpl
mpl.use("Agg")

import os
import sys
import numpy as np
import netCDF4
import operator
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

download_path = '/home/vmcd/' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

filebase = download_path + 'FYSP_clouds_archive/CAM4/'
outfileloc = download_path + 'temp_data/' # this is the location to save the processed netcdf files to

current = '0.775'
cc = '0775'

# SWCF fraction variable
field = 'SWCF'
outfilebase = 'c4_swcf_'
casenames = {'07','0725','075','0775', '08','0825','085','0875','09','0925','095','0975','10', '1025', '105', '1075','11'}

# 1.0 case
outfile_10 = outfileloc + outfilebase + '10.nc'
if not os.path.isfile(outfile_10):
	if os.path.isdir(outfileloc):
		infile = filebase +'cam4_10.nc'
		# calc cldlow global average per month
		syscall = r"//usr//bin//cdo timmean -seltimestep,21/40 -select,name="+field+" "+infile+ " " +outfile_10
		os.system(syscall)

for c in casenames:
	# calc swcf
	outfile_case = outfileloc+outfilebase+c+'.nc'
	# check directly if the file exists
	if not os.path.isfile(outfile_case):
		if os.path.isdir(outfileloc):
			infile = filebase +'cam4_' + c +'.nc'
			# calc cldlow global average per month
			syscall = r"//usr//bin//cdo timmean -seltimestep,21/40 -select,name="+field+" "+infile+ " " +outfile_case
			os.system(syscall)

control = outfile_10
if os.path.isfile(control):
	dsyear = netCDF4.Dataset(control)
	control_swcf = dsyear.variables[field][:]
	dsyear.close()


	#plot the data
	dsloc = outfileloc + outfilebase + cc +'.nc'
	if os.path.isfile(dsloc):
		# open the merged file and get out some variables
		dsyear = netCDF4.Dataset(dsloc)
		lons = dsyear.variables['lon'][:]
		lats = dsyear.variables['lat'][:]
		swcf = dsyear.variables[field][:]
		swcf_units = dsyear.variables[field].units
		dsyear.close() #close the file

		swcf_diff = list(map(operator.sub, swcf, control_swcf))

		#create plot
		fig = plt.figure()

		# setup the map
		m = Basemap(lat_0=0,lon_0=0)
		m.drawcoastlines()
		m.drawcountries()
	
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)

		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(swcf_diff), cmap='RdBu_r', latlon='True', vmin=-60, vmax=60, rasterized=True)

		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# Add Colorbar
		cbar = m.colorbar(cs, location='bottom', pad="10%")
		cbar.set_label(swcf_units)	
		
		plt.title('Shortwave Cloud Forcing: ' + r'$\mathsf{S/S_0}$'+' = '+ current)
		plt.show()

		fig.savefig('swcf_map_diff_'+cc+'.pdf', bbox_inches='tight')


