"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""
import os
import sys
import numpy as np
import netCDF4
from map_functions import get_map_field
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap

# cloud climatology
cloudfields= ['CLDHGH',
	'ICLDTWP',
	#'FICE',
	'LWCF',
	'CLDLOW',
	'ICLDTWP',
	#'FICE', 
	'SWCF']

cloudcmaps=['bone', 'BrBG', 'BuPu', 'BrBG', 'PuRd', 'RdBu_r', 'bone', 'BrBG', 'BuPu', 'BrBG', 'OrRd_r', 'RdBu_r']

cloudfilenames = ['map_annual_average',
	'tothighwp_map',
	#'icehighcloud',
	'map_annual_average',
	'map_annual_average',
	'totlowwp_map',
	#'icelowcloud',
	'map_annual_average']

cloudletters = ['a', 'b', 'c', 'd', 'e', 'f']
cloudheadings = ['High Cloud Fraction', 
	'Total High Cloud Water Path', 
	#'High CLoud Ice Fraction',
	'Longwave Cloud Forcing', 
	'Low Cloud Fraction', 
	'Total Low Cloud Water Path',
	#'Low Cloud Ice Fraction', 
	'Shortwave Cloud Forcing']

cloudaxislabels = [r'$\mathsf{Fraction}$', 
	r'$\mathsf{g/m^2}$', 
	#r'$\mathrm{Fraction}$', 
	r'$\mathsf{W/m^2}$', 
	r'$\mathsf{Fraction}$', 
	r'$\mathsf{g/m^2}$', 
	#r'$\mathrm{Fraction}$', 
	r'$\mathsf{W/m^2}$']

cloudvmins = [0,-0.5, 0, -300, 0, -60, 0, -0.5, 0,-300, -110, -60]
cloudvmaxs = [1, 0.5, 600, 300, 100, 60, 1, 0.5, 700, 300, 0, 60]



filebase = '/nazko/home/shared/vmcd/modelOutput/gcmexpt/CAM4/'

present = '1.0'
eight = '0.8'

#create plot
fig = plt.figure()

# container with 2 rows of 2 columns, first column is grid of absolute value plots, second column is diff plots. First row is cloud climatology, second row is model climatology
outer_grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.1, width_ratios=(2,1))

# first two columns, absolute value plots
cldabsgrid = gridspec.GridSpecFromSubplotSpec(6, 3, subplot_spec=outer_grid[0], wspace=0.0, hspace=0.2, width_ratios=(25,25,1))

# third colum, anomaly plots
clddiffgrid = gridspec.GridSpecFromSubplotSpec(6, 2, subplot_spec=outer_grid[1], wspace=0.0, hspace=0.2, width_ratios=(35,1))	



# -------------------------CLOUD CLIMATOLOGY -------------------------


# keep track of which field/row we're on
n=0
# keep track of which gridspace/column we're plotting in for abs val
a = 0
# keep track of which gridspace/column we're plotting in for diff		
d = 0
# keep track of which vmin/max we're on
v = 0
for p in cloudfields:
	f = cloudfilenames[n]
	cloudfield = cloudfields[n]
	presentcase = filebase + present+'/'+f+'.nc'
	eightcase = filebase + eight +'/'+f+'.nc'

	
	#plot the data - PRESENT
	#ax = grid[a]

	ax = fig.add_subplot(cldabsgrid[a])
	a=a+1

	if p == 'ICLDTWP':
		ds = netCDF4.Dataset(presentcase)
		presfld = ds.variables[p][:]
		lats = ds.variables['lat'][:]
		lons = ds.variables['lon'][:]

		ds.close() #close the file	

	else:
		output = get_map_field(p, present)

		presfld = output[0]
		lats = output[1]
		lons = output[2]


	# setup the map
	m = Basemap(lat_0=0,lon_0=0, ax=ax)
	m.drawcoastlines()
	m.drawcountries()
		
	# Create 2D lat/lon arrays for Basemap
	lon2d, lat2d = np.meshgrid(lons, lats)
	
	# Plot the data
	cs = m.pcolormesh(lon2d,lat2d,np.squeeze(presfld), cmap=cloudcmaps[v], latlon='True', vmin=cloudvmins[v], vmax=cloudvmaxs[v])
	
	# This is the fix for the white lines between contour levels
	cs.set_edgecolor("face")
		
	# add letter annotation
	plt.text(-0.10, 1.0, cloudletters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)

	# add heading
	plt.text(0.65, 1.05, cloudheadings[n], fontsize=12, transform=ax.transAxes)


	#plot the data - EIGHT
	ax = fig.add_subplot(cldabsgrid[a])
	a=a+1

	if p == 'ICLDTWP':
		ds = netCDF4.Dataset(eightcase)
		efld = ds.variables[p][:]
		lats = ds.variables['lat'][:]
		lons = ds.variables['lon'][:]

		ds.close() #close the file
	else:
		output = get_map_field(p, eight)

		efld = output[0]
		lats = output[1]
		lons = output[2]
	
	# setup the map
	m = Basemap(lat_0=0,lon_0=0, ax=ax)
	m.drawcoastlines()
	m.drawcountries()
		
	# Create 2D lat/lon arrays for Basemap
	lon2d, lat2d = np.meshgrid(lons, lats)

	# Plot 
	cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld), cmap=cloudcmaps[v], latlon='True', vmin=cloudvmins[v], vmax=cloudvmaxs[v])
	v = v+1
	
	# This is the fix for the white lines between contour levels
	cs.set_edgecolor("face")

	# plot the colorbar - ABS value
	ax = fig.add_subplot(cldabsgrid[a])
	a=a+1
	cb = plt.colorbar(cs, cax=ax)

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
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)
	
		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld)-np.squeeze(presfld), cmap=cloudcmaps[v], latlon='True', vmin=cloudvmins[v], vmax=cloudvmaxs[v])
		v = v+1
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# plot the colorbar - DIFF value
		ax = fig.add_subplot(clddiffgrid[d])
		d=d+1
		cb = plt.colorbar(cs, cax=ax, label=cloudaxislabels[n])

		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()

	# go to next field/row
	n=n+1


# -----------------------------

plt.show()

fig.savefig("figure4.pdf", bbox_inches='tight')

