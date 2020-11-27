#!/usr/bin/env python3
"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""
import matplotlib
matplotlib.use("Agg")

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap
import processing_functions as pf


# ------------------------------------------------------------------------
# change this section to match where you downloaded the model output files 
# ------------------------------------------------------------------------

download_path = '/home/vmcd/' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

filebase = download_path + 'FYSP_clouds_archive/CAM5/'
outfileloc = download_path + 'temp_data/' # this is the location to save the processed netcdf files to

# ------------------------------------

fields = 'LWCF,SWCF'
pf.cloud_forcing_all(filebase, outfileloc, fields, 'cam5')


outfilebase = 'c5_cloudforcing'

outer_grid = gridspec.GridSpec(2, 1, wspace=0.1, hspace=0.1)

#create plot
fig = plt.figure(figsize=(8.5, 9))
i=0
n=0
l=0

swcf = 'SWCF'
lwcf = 'LWCF'

present = '_10'
casenames = ['_105','_1025','_10','_0975','_095','_0925','_09' ]

while i<2:
	inner_grid = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=outer_grid[i], wspace=0.0, hspace=0.3)	
	if i==0:
		field = swcf
		outfilebase=outfilebase
	else:
		field = lwcf
		outfilebase=outfilebase
		n=0
		l=1
		

	for y in casenames:
		CASENAME = casenames[n]

		ax = fig.add_subplot(inner_grid[n])		
		
		n=n+1
	
		#plot the data
		presentcase = outfileloc + outfilebase + present +'.nc'
		currentcase = outfileloc + outfilebase + CASENAME +'.nc'

		ds = netCDF4.Dataset(presentcase)
		presfld = ds.variables[field][:]
		lats = ds.variables['lat'][:]
		lons = ds.variables['lon'][:]

		ds = netCDF4.Dataset(currentcase)
		cfld = ds.variables[field][:]
		lats = ds.variables['lat'][:]
		lons = ds.variables['lon'][:]

	
		# setup the map
		m = Basemap(lat_0=0,lon_0=0,ax=ax)
		m.drawcoastlines()
		m.drawcountries()
		parallels = [-45, 0, 45]
		meridians = [-90., 0., 90.]

		if n==1 or n==4 or n==7:
			m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
		else:
			m.drawparallels(parallels, labels=[False ,False,False, False], fontsize=6)

		m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
	
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)
	
		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(cfld) - np.squeeze(presfld), cmap='RdBu_r', latlon='True', vmin=-60, vmax=60, rasterized=True)
		fig.add_subplot(ax)
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")
		
		# add letter annotation
		if l==0:
			plt.text(-0.10, 1.0, "a", fontsize=10, fontweight="bold", transform=ax.transAxes)
		elif l==1:
			plt.text(-0.10, 1.0, "b", fontsize=10, fontweight="bold", transform=ax.transAxes)

		l=l+2
	
		ax.set_title(CASENAME[1]+'.'+CASENAME[-1], fontsize=10)
	i=i+1
	
# Add Colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cb = fig.colorbar(cs, cax=cbar_ax)
cb.set_label(label=r'$\mathsf{W/m^2}$', fontsize=6)

cb.ax.tick_params(labelsize=6) 
cb.update_ticks()

plt.show()

fig.savefig("ED_figure6.pdf", bbox_inches='tight')

