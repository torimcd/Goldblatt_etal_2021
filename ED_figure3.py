#!/usr/bin/env python3

"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

This script creates maps and anomaly maps of the model cloud climatology in CAM5.

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

download_path = '' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

filebase = download_path + 'FYSP_clouds_archive/CAM5/'
outfileloc = download_path + 'temp_data/' # this is the location to save the processed netcdf files to

# the fields we want to average for our plots  - these must not depend on pressure level
fields = 'CLDHGH,CLDLOW,LHFLX,LWCF,PRECT,SHFLX,SWCF,TS'

# process the fields we're plotting
pf.map_annual_average(filebase, outfileloc, 'cam5', fields) # averages fields over years 31-60, retaining location so can be plotted in map view


# cloud climatology
cloudfields= ['CLDHGH', 'LWCF', 'CLDLOW', 'SWCF']

cloudcmaps=['bone', 'BrBG', 'PuRd', 'RdBu_r', 'bone', 'BrBG', 'OrRd_r', 'RdBu_r']

cloudfilenames = ['c5_map_annual_average', 'c5_map_annual_average', 'c5_map_annual_average', 'c5_map_annual_average']

cloudletters = ['a', 'b', 'c', 'd']
cloudheadings = ['High Cloud Fraction',	'Longwave Cloud Forcing', 'Low Cloud Fraction', 'Shortwave Cloud Forcing']

cloudaxislabels = [r'$\mathrm{Fraction}$', r'$\mathrm{W/m^2}$', r'$\mathrm{Fraction}$', r'$\mathrm{W/m^2}$']

cloudvmins = [0,-0.5, 0, -60, 0, -0.5, -110, -60]
cloudvmaxs = [1, 0.5, 100, 60, 1, 0.5, 0, 60]


#create figure - use figsize=(8.5, 9.5) to make bigger
#fig = plt.figure(figsize=(3.46457, 4.48356))
fig = plt.figure(figsize=(7.08661, 5.2107))


# container with 2 rows of 2 columns, first column is grid of absolute value plots, second column is diff plots. First row is cloud climatology, second row is model climatology
outer_grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.1, width_ratios=(2,1))

# first two columns, absolute value plots
cldabsgrid = gridspec.GridSpecFromSubplotSpec(4, 3, subplot_spec=outer_grid[0], wspace=0.0, hspace=0.45, width_ratios=(15,15,1))

# third colum, anomaly plots
clddiffgrid = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=outer_grid[1], wspace=0.0, hspace=0.45, width_ratios=(25,1))	



# -------------------------CLOUD CLIMATOLOGY -------------------------


# keep track of which field/row we're on
n=0
# keep track of which gridspace/column we're plotting in for abs val
a = 0
# keep track of which gridspace/column we're plotting in for diff		
d = 0
# keep track of which vmin/max we're on
v = 0

present = '_10'
eight = '_09'

for p in cloudfields:
	f = cloudfilenames[n]
	cloudfield = cloudfields[n]

	presentcase = outfileloc + f + present +'.nc'
	eightcase = outfileloc + f + eight +'.nc'
	
	# plot the data - PRESENT
	ax = fig.add_subplot(cldabsgrid[a])
	a=a+1

	ds = netCDF4.Dataset(presentcase)
	lons = ds.variables['lon'][:]
	lats = ds.variables['lat'][:]
	presfld = ds.variables[cloudfield][:]
	units = ds.variables[cloudfield].units
	
	ds.close() #close the file

	# setup the map
	m = Basemap(lat_0=0,lon_0=0, ax=ax)
	m.drawcoastlines()
	m.drawcountries()
	parallels = [-45, 0, 45]
	meridians = [-90., 0., 90.]
	m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
	m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
		
	# Create 2D lat/lon arrays for Basemap
	lon2d, lat2d = np.meshgrid(lons, lats)
	
	# Plot the data
	cs = m.pcolormesh(lon2d,lat2d,np.squeeze(presfld), cmap=cloudcmaps[v], latlon='True', vmin=cloudvmins[v], vmax=cloudvmaxs[v], rasterized=True)
	
	# This is the fix for the white lines between contour levels
	cs.set_edgecolor("face")
		
	# add letter annotation
	plt.text(-0.10, 1.0, cloudletters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)

	# add heading
	plt.text(0.65, 1.05, cloudheadings[n], fontsize=7, transform=ax.transAxes)


	#plot the data - EIGHT
	ax = fig.add_subplot(cldabsgrid[a])
	a=a+1

	ds = netCDF4.Dataset(eightcase)
	lons = ds.variables['lon'][:]
	lats = ds.variables['lat'][:]
	efld = ds.variables[cloudfield][:]
	units = ds.variables[cloudfield].units
	
	ds.close() #close the file
	
	# setup the map
	m = Basemap(lat_0=0,lon_0=0, ax=ax)
	m.drawcoastlines()
	m.drawcountries()
	parallels = [-45, 0, 45]
	meridians = [-90., 0., 90.]
	m.drawparallels(parallels, labels=[False ,False,False, False], fontsize=6)
	m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
		
	# Create 2D lat/lon arrays for Basemap
	lon2d, lat2d = np.meshgrid(lons, lats)

	# Plot 
	cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld), cmap=cloudcmaps[v], latlon='True', vmin=cloudvmins[v], vmax=cloudvmaxs[v], rasterized=True)
	v = v+1
	
	# This is the fix for the white lines between contour levels
	cs.set_edgecolor("face")

	# plot the colorbar - ABS value
	ax = fig.add_subplot(cldabsgrid[a])
	a=a+1
	cb = plt.colorbar(cs, cax=ax)

	cb.ax.tick_params(labelsize=6) 
	tick_locator = ticker.MaxNLocator(nbins=5)
	cb.locator = tick_locator
	cb.update_ticks()

	#plot the data - DIFF
	ax = fig.add_subplot(clddiffgrid[d])
	d=d+1
	if os.path.isfile(eightcase):
	
		# setup the map
		m = Basemap(lat_0=0,lon_0=0, ax=ax)
		m.drawcoastlines()
		m.drawcountries()
		parallels = [-45, 0, 45]
		meridians = [-90., 0., 90.]
		m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
		m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)
	
		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld)-np.squeeze(presfld), cmap=cloudcmaps[v], latlon='True', vmin=cloudvmins[v], vmax=cloudvmaxs[v], rasterized=True)
		v = v+1
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# plot the colorbar - DIFF value
		ax = fig.add_subplot(clddiffgrid[d])
		d=d+1
		cb = plt.colorbar(cs, cax=ax)

		cb.set_label(label=cloudaxislabels[n], fontsize=6)
		cb.ax.tick_params(labelsize=6) 
		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()

	# go to next field/row
	n=n+1


# -----------------------------

plt.show()

fig.savefig("figures_ED/ED_figure3.pdf", format='pdf', bbox_inches='tight')

