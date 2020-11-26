#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap


# model climatology
climfields= ['TS', 'LHFLX', 'SHFLX']

climcmaps=['plasma', 'RdBu_r', 'YlOrBr', 'RdBu_r', 'YlOrBr', 'RdBu_r']

climfilenames = ['map_annual_average', 'map_annual_average', 'map_annual_average', 'mapvertvelc5']
climletters = ['a', 'b','c','d','e','f','g']
climheadings = ['Surface Temperature', 'Latent Heat Flux', 'Sensible Heat Flux', 'Total Precipitation Rate', 'Vertical Velocity at 700 hPa',  'Lower Tropospheric Stability', 'Estimated Inversion Strength']


climvmins = [220,-10, 0, -75, 0, -75, -12, -6]
climvmaxs = [320, 10, 250, 75, 250, 75, 12, 6]

eis_lcl_file = 'eis_map_lcl.nc'
eis_qs850_file = 'eis_map_qs850.nc'
eis_ = ''

filebase = '/nazko/home/shared/vmcd/modelOutput/gcmexpt/CAM5/'

present = '1.0'
eight = '0.9'

#create plot
fig = plt.figure(figsize=(8.5, 11))

# container with 2 rows of 2 columns, first column is grid of absolute value plots, second column is diff plots. First row is cloud climatology, second row is model climatology
outer_grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.1, width_ratios=(2,1))


# first two columns, absolute value plots
climabsgrid = gridspec.GridSpecFromSubplotSpec(7, 3, subplot_spec=outer_grid[0], wspace=0.0, hspace=0.4, width_ratios=(25,25,1))

# third colum, anomaly plots
climdiffgrid = gridspec.GridSpecFromSubplotSpec(7, 2, subplot_spec=outer_grid[1], wspace=0.0, hspace=0.4, width_ratios=(35,1))	



# ----------------------- MODEL CLIMATOLOGY -----------------------------

# keep track of which field/row we're on
n=0
# keep track of which gridspace/column we're plotting in for abs val
a = 0
# keep track of which gridspace/column we're plotting in for diff		
d = 0
# keep track of which vmin/max we're on
v = 0



for p in climfields:
	f = climfilenames[n]

	climfield = climfields[n]
	presentcase = filebase + present+'/'+f+'.nc'
	eightcase = filebase + eight+'/'+f+'.nc'

	if f == 'mapvertvelc5':
		presentcase = '/home/vmcd/projects/output/mapvertvelc5_1.0.nc'
		eightcase = '/home/vmcd/projects/output/mapvertvelc5_0.9.nc'


	#plot the data - PRESENT
	#ax = grid[a]

	ax = fig.add_subplot(climabsgrid[a])
	a=a+1

	if os.path.isfile(presentcase):
		# open the file and get out some variables
		dsyear = netCDF4.Dataset(presentcase)
		lons = dsyear.variables['lon'][:]
		lats = dsyear.variables['lat'][:]
		presfld = dsyear.variables[climfield][:]
		units = dsyear.variables[climfield].units
	
		dsyear.close() #close the file

		if units == 'Pa/s':
			presfld = presfld*100
			units = 'hPa/s'

		# setup the map
		m = Basemap(lat_0=0,lon_0=0, ax=ax)
		m.drawcoastlines()
		m.drawcountries()
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)
	
		# Plot the data
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(presfld), cmap=climcmaps[v], latlon='True', vmin=climvmins[v], vmax=climvmaxs[v], rasterized=True)
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")
		
		# add letter annotation
		plt.text(-0.10, 1.0, climletters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)
		# add heading
		plt.text(0.65, 1.05, climheadings[n], fontsize=12, transform=ax.transAxes)

	#plot the data - EIGHT
	ax = fig.add_subplot(climabsgrid[a])
	a=a+1

	if os.path.isfile(eightcase):
		# open the file and get out some variables
		dsyear = netCDF4.Dataset(eightcase)
		lons = dsyear.variables['lon'][:]
		lats = dsyear.variables['lat'][:]
		efld = dsyear.variables[climfield][:]
		units = dsyear.variables[climfield].units
	
	
		if units == 'Pa/s':
			efld = efld*100
			units = 'hPa/s'

		# setup the map
		m = Basemap(lat_0=0,lon_0=0, ax=ax)
		m.drawcoastlines()
		m.drawcountries()
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)

		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld), cmap=climcmaps[v], latlon='True', vmin=climvmins[v], vmax=climvmaxs[v], rasterized=True)
		v = v+1
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# plot the colorbar - ABS value
		ax = fig.add_subplot(climabsgrid[a])
		a=a+1
		cb = plt.colorbar(cs, cax=ax)

		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()
		

	#plot the data - DIFF
	ax = fig.add_subplot(climdiffgrid[d])
	d=d+1
	if os.path.isfile(eightcase):
	
		# setup the map
		m = Basemap(lat_0=0,lon_0=0, ax=ax)
		m.drawcoastlines()
		m.drawcountries()
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)
	
		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld)-np.squeeze(presfld), cmap=climcmaps[v], latlon='True', vmin=climvmins[v], vmax=climvmaxs[v], rasterized=True)
		v = v+1
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# plot the colorbar - DIFF value
		ax = fig.add_subplot(climdiffgrid[d])
		d=d+1
		cb = plt.colorbar(cs, cax=ax, label=units)

		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()

	# go to next field/row
	n=n+1


# Special Case PREC -----------------------

presentcase = filebase + present+'/map_annual_average.nc'
eightcase = filebase + eight+'/map_annual_average.nc'


prectotalpres = []
prectotal = []
prectotaldiff = []

if os.path.isfile(presentcase):
	# open the file and get out some variables
	dsyear = netCDF4.Dataset(presentcase)
	lons = dsyear.variables['lon'][:]
	lats = dsyear.variables['lat'][:]
	present_precc = dsyear.variables['PRECC'][:]
	present_precl = dsyear.variables['PRECL'][:]

	prectotalpres = present_precc + present_precl

if os.path.isfile(eightcase):
	# open the file and get out some variables
	dsyear = netCDF4.Dataset(eightcase)
	lons = dsyear.variables['lon'][:]
	lats = dsyear.variables['lat'][:]
	eight_precc = dsyear.variables['PRECC'][:]
	eight_precl = dsyear.variables['PRECL'][:]

	prectotal = eight_precc + eight_precl


prectotaldiff = prectotal - prectotalpres

##    -----------------------------


## Special Case Omega -----------

presentcase = '/home/vmcd/projects/output/mapvertvelc5_1.0.nc'
eightcase = '/home/vmcd/projects/output/mapvertvelc5_0.9.nc'


if os.path.isfile(presentcase):
	# open the file and get out some variables
	dsyear = netCDF4.Dataset(presentcase)
	lons = dsyear.variables['lon'][:]
	lats = dsyear.variables['lat'][:]
	pres_omega = dsyear.variables['OMEGA'][:]
	units = dsyear.variables['OMEGA'].units
	
	dsyear.close() #close the file

	if units == 'Pa/s':
		pres_omega = pres_omega*100
		units = 'hPa/s'

if os.path.isfile(eightcase):
	# open the file and get out some variables
	dsyear = netCDF4.Dataset(eightcase)
	lons = dsyear.variables['lon'][:]
	lats = dsyear.variables['lat'][:]
	e_omega = dsyear.variables['OMEGA'][:]
	units = dsyear.variables['OMEGA'].units
	
	
	if units == 'Pa/s':
		e_omega = e_omega*100
		units = 'hPa/s'
	
omega_diff = e_omega - pres_omega

##    ----------------------------


# Special Case - EIS ------------
outfilebase = 'eis_map_'
lts_outfilebase = 'lts_map_'


qs850_p = filebase + present+'/'+outfilebase+'_qs850.nc'
temp700_p = filebase + present+'/'+outfilebase+'_temp700.nc'
tempsurf_p = filebase + present+'/'+outfilebase+'_tempsurf.nc'
tempsum_p = filebase + present+'/'+outfilebase+'_tempsum.nc'
z700_p = filebase + present+'/'+outfilebase+'_z700.nc'
lcl_p = filebase + present+'/'+outfilebase+'_lcl.nc'

qs850_e = filebase + eight+'/'+outfilebase+'_qs850.nc'
temp700_e = filebase + eight+'/'+outfilebase+'_temp700.nc'
tempsurf_e = filebase + eight+'/'+outfilebase+'_tempsurf.nc'
tempsum_e = filebase + eight+'/'+outfilebase+'_tempsum.nc'
z700_e = filebase + eight+'/'+outfilebase+'_z700.nc'

lcl_e = filebase + eight+'/'+outfilebase+'_lcl.nc'

lcl = []
z700 = []
lts = []
qs850 = []
tempsum = []
pal_e = []
eis_e = []
lons = []
lats = []

lcl_pres = []
z700_pres = []
lts_pres = []
qs850_pres = []
tempsum_pres = []
pal_pres = []
eis_pres = []
lons_pres = []
lats_pres = []

	
# GET LCL eight
if os.path.isfile(lcl_e):
	# open the file and get out the variable
	ds = netCDF4.Dataset(lcl_e)
	lcl = ds.variables['lcl'][:]
	ds.close() #close the file

# GET LCL Present
if os.path.isfile(lcl_p):
	# open the file and get out the variable
	ds = netCDF4.Dataset(lcl_p)
	lcl_pres = ds.variables['lcl'][:]
	ds.close() #close the file





#GET z700 eight
if os.path.isfile(z700_e):
	# open the file and get out the variable
	ds = netCDF4.Dataset(z700_e)
	z700 = ds.variables['Z3'][:]
	ds.close() #close the file

#GET z700 Present
if os.path.isfile(z700_p):
	# open the file and get out the variable
	ds = netCDF4.Dataset(z700_p)
	z700_pres = ds.variables['Z3'][:]
	ds.close() #close the file




#GET lts eight
dsloc = filebase + eight+'/'+lts_outfilebase+'sub.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	lts = ds.variables['lts'][:]
	ds.close() #close the file

#GET lts Present
dsloc = filebase + present+'/'+lts_outfilebase+'sub.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	lts_pres = ds.variables['lts'][:]
	ds.close() #close the file




#GET qs850 eight
if os.path.isfile(qs850_e):
	# open the file and get out the variable
	ds = netCDF4.Dataset(qs850_e)
	qs850 = ds.variables['smr'][:]
	ds.close() #close the file

#GET qs850 Present
if os.path.isfile(qs850_p):

	# open the file and get out the variable
	ds = netCDF4.Dataset(qs850_p)
	qs850_pres = ds.variables['smr'][:]
	ds.close() #close the file	





#GET tempsum eight
if os.path.isfile(tempsum_e):

	# open the file and get out the variable
	ds = netCDF4.Dataset(tempsum_e)
	tempsum = ds.variables['T'][:]
	lons = ds.variables['lon'][:]
	lats = ds.variables['lat'][:]
	ds.close() #close the file

#GET tempsum Present
if os.path.isfile(tempsum_p):

	# open the file and get out the variable
	ds = netCDF4.Dataset(tempsum_p)
	tempsum_pres = ds.variables['T'][:]
	lons_pres = ds.variables['lon'][:]
	lats_pres = ds.variables['lat'][:]
	ds.close() #close the file



pal_e = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	
eis_e = lts - pal_e*((z700/1000)-lcl)

	
pal_present = (9.9)*(1-(1+2450000*qs850_pres/(287.058*(tempsum_pres/2)))/(1+(2450000**2)*qs850_pres/(993*461.4*((tempsum_pres/2)**2))))	
eis_present = lts_pres - pal_present*((z700_pres/1000)-lcl_pres)

ltsdiff = lts - lts_pres
eis_diff = eis_e - eis_present







# Plot PREC Present  **************************
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(prectotalpres), cmap='PuBuGn', latlon='True', vmin=0, vmax=0.00000015, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=12, transform=ax.transAxes)

# Plot PREC Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(prectotal), cmap='PuBuGn', latlon='True', vmin=0, vmax=0.00000015, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


# Plot PREC Diff **************************** 
ax = fig.add_subplot(climdiffgrid[d])
d=d+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(prectotaldiff), cmap='BrBG', latlon='True', vmin=-0.00000005,vmax=0.00000005, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax, label='K')

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

n = n+1

# -----------------------------



# Plot OMEGA Present  **************************
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(pres_omega), cmap='RdGy', latlon='True', vmin=-10, vmax=10, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=12, transform=ax.transAxes)



# Plot OMEEGA Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(e_omega), cmap='RdGy', latlon='True', vmin=-10, vmax=10, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


# Plot OMEGA Diff **************************** 
ax = fig.add_subplot(climdiffgrid[d])
d=d+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(omega_diff), cmap='RdGy', latlon='True', vmin=-4,vmax=4, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax, label='hPa/s')

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

n = n+1

# -----------------------------



# Plot LTS Present  **************************
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(lts_pres), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=12, transform=ax.transAxes)

# Plot LTS Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(lts), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


# Plot LTS Diff **************************** 
ax = fig.add_subplot(climdiffgrid[d])
d=d+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(ltsdiff), cmap='PiYG_r', latlon='True', vmin=-15,vmax=15, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax, label='K')

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

n = n+1

# -----------------------------

# Plot EIS Present  **************************
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(eis_present), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=12, transform=ax.transAxes)

# Plot EIS Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(eis_e), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


# Plot EIS Diff **************************** 
ax = fig.add_subplot(climdiffgrid[d])
d=d+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(eis_diff), cmap='PiYG_r', latlon='True', vmin=-15,vmax=15, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax, label='K')

tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()



# -----------------------------

plt.show()

fig.savefig("figure2CAM5.pdf", bbox_inches='tight')

