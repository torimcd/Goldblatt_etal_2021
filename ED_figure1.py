#!/usr/bin/env python3

"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

This script creates maps and anomaly maps of the model climatology for CAM5 model output.

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
fields = 'CLDHGH,CLDLOW,LHFLX,LWCF,PRECC,PRECL,SHFLX,SWCF,TS'

# process the fields we're plotting
pf.map_annual_average(filebase, outfileloc, 'cam5', fields) # averages fields over years 31-60, retaining location so can be plotted in map view
pf.map_vert_velocity(filebase, outfileloc, 'cam5') # the same as above but for vertical velocity selected at 700 hPa
pf.map_prep_lts(filebase, outfileloc, 'cam5') # extracts and calculates LTS and saves to new file
pf.map_prep_eis(filebase, outfileloc, 'cam5') # extracts and calculates EIS and saves to new file


# model climatology
climfields= ['TS', 'LHFLX', 'SHFLX']

climcmaps=['plasma', 'RdBu_r', 'YlOrBr', 'RdBu_r', 'YlOrBr', 'RdBu_r']

climfilenames = ['c5_map_annual_average', 'c5_map_annual_average', 'c5_map_annual_average']
climletters = ['a', 'b','c','d','e','f','g']
climheadings = ['Surface Temperature', 'Latent Heat Flux', 'Sensible Heat Flux', 'Total Precipitation Rate', 'Vertical Velocity at 700 hPa',  'Lower Tropospheric Stability', 'Estimated Inversion Strength']

# set the labels on the colorbars
units_all = [r'$\mathsf{K}$', r'$\mathsf{W/m^2}$', r'$\mathsf{W/m^2}$', r'$\mathsf{m/yr}$', r'$\mathsf{hPa/s}$', r'$\mathsf{K}$', r'$\mathsf{K}$']


climvmins = [220,-10, 0, -75, 0, -75, -12, -6]
climvmaxs = [320, 10, 250, 75, 250, 75, 12, 6]

# calculating EIS is more involved so there are separate processed files to read for these maps
eis_lcl_file = 'c5_eis_map_lcl.nc'
eis_qs850_file = 'c5_eis_map_qs850.nc'
eis_ = ''

#create figure 
fig = plt.figure(figsize=(7.08661, 9.17))

# set up container 
outer_grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.1, width_ratios=(2,1))

# first two columns, absolute value plots
climabsgrid = gridspec.GridSpecFromSubplotSpec(7, 3, subplot_spec=outer_grid[0], wspace=0.0, hspace=0.4, width_ratios=(15,15,1))

# third colum, anomaly plots
climdiffgrid = gridspec.GridSpecFromSubplotSpec(7, 2, subplot_spec=outer_grid[1], wspace=0.0, hspace=0.4, width_ratios=(25,1))	



# ----------------------------
# Model Climatology
#-----------------------------

# keep track of which field/row we're on
n=0
# keep track of which gridspace/column we're plotting in for abs val
a = 0
# keep track of which gridspace/column we're plotting in for diff		
d = 0
# keep track of which vmin/max we're on for the colorbar
v = 0

present = '_10'
eight = '_09' # lowest S/S0 CAM5 run is 0.9

# get the data for the first five fields
for p in climfields:
	f = climfilenames[n]
	climfield = climfields[n]

	# get out the data for the 1.0 S/So and 0.8 S/So
	presentcase = outfileloc + f + present +'.nc'
	eightcase = outfileloc + f + eight +'.nc'

	
	# add the subplot for this field
	ax = fig.add_subplot(climabsgrid[a])
	a=a+1

	# plot the S/So = 1.0 case
	if os.path.isfile(presentcase):
		# open the file and get out some variables
		dsyear = netCDF4.Dataset(presentcase)
		lons = dsyear.variables['lon'][:]
		lats = dsyear.variables['lat'][:]
		presfld = dsyear.variables[climfield][:] # presfld is the present day 1.0 field
		units = dsyear.variables[climfield].units
	
		dsyear.close() #close the file

		# scale it so the units consistent
		if units == 'm/s':
			presfld = presfld*31540000
			units = 'm/yr'

		if units == 'Pa/s':
			presfld = presfld*100
			units = 'hPa/s'

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
	
		# Plot the data -- rasterized=true makes the image much smaller so it renders quickly
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(presfld), cmap=climcmaps[v], latlon='True', vmin=climvmins[v], vmax=climvmaxs[v], rasterized=True)
	
		# This removes white lines between contour levels
		cs.set_edgecolor("face")
		
		# add letter annotation
		plt.text(-0.10, 1.0, climletters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)
		# add heading
		plt.text(0.65, 1.05, climheadings[n], fontsize=7, transform=ax.transAxes)

	# plot the S/So = 0.8 case
	ax = fig.add_subplot(climabsgrid[a])
	a=a+1

	if os.path.isfile(eightcase):
		# open the file and get out some variables
		dsyear = netCDF4.Dataset(eightcase)
		lons = dsyear.variables['lon'][:]
		lats = dsyear.variables['lat'][:]
		efld = dsyear.variables[climfield][:]
		units = dsyear.variables[climfield].units

		if units == 'm/s':
			efld = efld*31540000
			units= 'm/yr'

		if units == 'Pa/s':
			efld = efld*100
			units = 'hPa/s'
	
	
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
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld), cmap=climcmaps[v], latlon='True', vmin=climvmins[v], vmax=climvmaxs[v], rasterized=True)
		v = v+1
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# plot the colorbar - ABS value
		ax = fig.add_subplot(climabsgrid[a])
		a=a+1
		cb = plt.colorbar(cs, cax=ax)

		cb.ax.tick_params(labelsize=6) 
		cb.ax.yaxis.offsetText.set(size=6)
		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()

	# plot the difference field, 0.8 - 1.0
	ax = fig.add_subplot(climdiffgrid[d])
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
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(efld)-np.squeeze(presfld), cmap=climcmaps[v], latlon='True', vmin=climvmins[v], vmax=climvmaxs[v], rasterized=True)
		v = v+1
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")

		# plot the colorbar - DIFF value
		ax = fig.add_subplot(climdiffgrid[d])
		d=d+1
		cb = plt.colorbar(cs, cax=ax)

		cb.set_label(label=units_all[n], fontsize=6)
		cb.ax.tick_params(labelsize=6) 
		cb.ax.yaxis.offsetText.set(size=6)
		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()

	# go to next field/row
	n=n+1



# ----------------------------
# Precip is a special case
#-----------------------------
presentcase = outfileloc + '/c5_map_annual_average' + present + '.nc'
eightcase = outfileloc +'/c5_map_annual_average' + eight + '.nc'


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


# ----------------------------
# Omega is a special case
#-----------------------------

presentcase = outfileloc + 'c5_map_vert_velocity' + present + '.nc'

eightcase = outfileloc + 'c5_map_vert_velocity' + eight + '.nc'


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



# ----------------------------
# EIS is a special case
#-----------------------------
outfilebase = 'c5_eis_map'
lts_outfilebase = 'c5_lts_map'


# all the variables we need were processed into separate files
qs850_p = outfileloc + outfilebase+ present + '_qs850.nc'
temp700_p = outfileloc + outfilebase+ present + '_temp700.nc'
tempsurf_p = outfileloc + outfilebase+ present + '_tempsurf.nc'
tempsum_p = outfileloc + outfilebase+ present + '_tempsum.nc'
z700_p = outfileloc + outfilebase+ present + '_z700.nc'
lcl_p = outfileloc + outfilebase+ present + '_lcl.nc'


qs850_e = outfileloc + outfilebase+ eight + '_qs850.nc'
temp700_e = outfileloc + outfilebase+ eight + '_temp700.nc'
tempsurf_e = outfileloc + outfilebase+ eight + '_tempsurf.nc'
tempsum_e = outfileloc + outfilebase+ eight + '_tempsum.nc'
z700_e = outfileloc + outfilebase+ eight + '_z700.nc'

lcl_e = outfileloc + outfilebase+ eight + '_lcl.nc'

# set up the arrays
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
dsloc = outfileloc + lts_outfilebase+ eight +'_lts.nc'
if os.path.isfile(dsloc):

	# open the file and get out the variable
	ds = netCDF4.Dataset(dsloc)
	lts = ds.variables['lts'][:]
	ds.close() #close the file

#GET lts Present
dsloc = outfileloc + lts_outfilebase+ present +'_lts.nc'
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



# ----------------------------
# Plot Precip
#-----------------------------
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(prectotalpres), cmap='PuBuGn', latlon='True', vmin=0, vmax=0.00000015, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=7, transform=ax.transAxes)

# Plot PREC Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(prectotal), cmap='PuBuGn', latlon='True', vmin=0, vmax=0.00000015, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
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
parallels = [-45, 0, 45]
meridians = [-90., 0., 90.]
m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(prectotaldiff), cmap='BrBG', latlon='True', vmin=-0.00000005,vmax=0.00000005, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax)

cb.set_label(label='K', fontsize=6)
cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

n = n+1


# ----------------------------
# Plot OMEGA
#-----------------------------
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(pres_omega), cmap='RdGy', latlon='True', vmin=-10, vmax=10, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=7, transform=ax.transAxes)



# Plot OMEEGA Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(e_omega), cmap='RdGy', latlon='True', vmin=-10, vmax=10, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
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
parallels = [-45, 0, 45]
meridians = [-90., 0., 90.]
m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(omega_diff), cmap='RdGy', latlon='True', vmin=-4,vmax=4, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax)

cb.set_label(label='hPa/s', fontsize=6)
cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

n = n+1




# ----------------------------
# Plot LTS
#-----------------------------
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(lts_pres), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=7, transform=ax.transAxes)

# Plot LTS Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(lts), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - ABS value
ax = fig.add_subplot(climabsgrid[a])
a=a+1
cb = plt.colorbar(cs, cax=ax)

cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
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
parallels = [-45, 0, 45]
meridians = [-90., 0., 90.]
m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
# Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lons, lats)
	
# Plot the data
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(ltsdiff), cmap='PiYG_r', latlon='True', vmin=-15,vmax=15, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")

# plot the colorbar - DIFF value
ax = fig.add_subplot(climdiffgrid[d])
d=d+1
cb = plt.colorbar(cs, cax=ax)

cb.set_label(label='K', fontsize=6)
cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

n = n+1



# ----------------------------
# PLOT EIS
#-----------------------------
ax = fig.add_subplot(climabsgrid[a])
a=a+1

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
cs = m.pcolormesh(lon2d,lat2d,np.squeeze(eis_present), cmap='PiYG_r', latlon='True', vmin=-40,vmax=40, rasterized=True)
	
# This is the fix for the white lines between contour levels
cs.set_edgecolor("face")
		
# add letter annotation
plt.text(-0.10, 1.0, climletters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)
# add heading
plt.text(0.65, 1.05, climheadings[n], fontsize=7, transform=ax.transAxes)

# Plot EIS Eight **************************** 
ax = fig.add_subplot(climabsgrid[a])
a=a+1

# setup the map
m = Basemap(lat_0=0,lon_0=0, ax=ax)
m.drawcoastlines()
m.drawcountries()
parallels = [-45, 0, 45]
meridians = [-90., 0., 90.]
m.drawparallels(parallels, labels=[False ,False,False, False], fontsize=5)
m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=5)
		
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

cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
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
parallels = [-45, 0, 45]
meridians = [-90., 0., 90.]
m.drawparallels(parallels, labels=[True ,False,False, False], fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,True], fontsize=6)
		
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

cb.ax.tick_params(labelsize=6) 
cb.ax.yaxis.offsetText.set(size=6)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

# -----------------------------


plt.show()

fig.savefig("ED_figure1.pdf", format='pdf', bbox_inches='tight')

