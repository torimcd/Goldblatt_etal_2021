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
import operator
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap
import processing_functions as pf

# enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'
download_path = '/Users/victoria/projects/data/fysp_clouds/' 

filebase = download_path + '/FYSP_clouds_archive/CAM4/'
filebase = download_path # FOR TESTING, REMOVE FOR DATA ARCHIVE

# process the fields we're plotting
pf.zonal_average(filebase, 'cam4') # averages fields zonally over years 31-60, retaining location so can be plotted in map view
pf.zonal_wbpt(filebase, 'cam4') # the same as above but for wetbulb potential temperature selected at 700 hPa


# variables
fields= ['T', 	# temperature
	'sigmaw', 	# wet bulb potential temperature
	'Q', 		# humidity
	'OMEGA']	# vertical velocity

# color maps
cmaps=['cool', 'cool','GnBu', 'RdGy']
cmaps_d=['RdBu_r','RdBu_r','BrBG','RdGy']

filenames = ['zonal_average_',
	'wbpt_',
	'zonal_average_',
	'zonal_average_']

letters = ['a', 'b', 'c', 'd', 'e', 'f']

vmins = [200, -50, 240, -10, 0, -30, -5, -1.8] 
vmaxs = [602, 52, 320, 11, 20, 31, 6, 1.8]

present = '1.0'
eight = '0.8'

casenames = {'1.1': '#fcae91','1.075': '#fb6a4a','1.05': '#cb181d','1.025': '#dd3497','1.0': 'black','0.975':'#8d2689','0.95': '#9366af','0.925': '#9061ee','0.9': '#4c4cff','0.875': '#0073e6','0.85': '#2c7fb8','0.825': '#41b6c4','0.8': '#97dde5','0.775':'#9fbfc8','0.75': '#aec3cb','0.725': '#c7cccd','0.7': '#90999c'}

#create plot
fig = plt.figure(figsize=(18, 12))

# container
outer_grid = gridspec.GridSpec(2, 1, wspace=0, hspace=0.2, height_ratios=(1,3))

midgrid_top = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=outer_grid[0], wspace=0.2, hspace=0.1)
midgrid_bottom = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=outer_grid[1], wspace=0.2, hspace=0.3,width_ratios=(2,1))

# first two columns, absolute value plots
absgrid = gridspec.GridSpecFromSubplotSpec(4, 3, subplot_spec=midgrid_bottom[0], wspace=0.05, hspace=0.1, width_ratios=(35,35,1))

# third colum, anomaly plots
diffgrid = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=midgrid_bottom[1], wspace=0.05, hspace=0.1, width_ratios=(35,1))	

n=0

ax1 = fig.add_subplot(midgrid_top[0])
ax2 = fig.add_subplot(midgrid_top[1])
#ax3 = fig.add_subplot(midgrid_top[1,:])

# Top row - zonal surface temp and diff -----------------------------------------------
for CASENAME in casenames.keys():
	dsloc_a = filebase+CASENAME+'/tszonmean.nc'
	if os.path.isfile(dsloc_a):

		# open the merged file and get out the variables
		ds = netCDF4.Dataset(dsloc_a)
		ts_a = ds.variables['TS'][:]
		lat = ds.variables['lat'][:]
		ds.close() #close the file
		ts_a = ts_a.flatten()

	dsloc_d = filebase+CASENAME+'/tszonmean_diff.nc'
	if os.path.isfile(dsloc_d):

		# open the merged file and get out the variables
		ds = netCDF4.Dataset(dsloc_d)
		ts_d = ds.variables['TS'][:]
		lat = ds.variables['lat'][:]
		ds.close() #close the file
		ts_d = ts_d.flatten()

	#plot the data
	ax1.plot(lat, ts_a, color=casenames[CASENAME], linewidth=1.5, label=CASENAME, rasterized=True)
	ax2.plot(lat, ts_d, color=casenames[CASENAME], linewidth=1.5, label=CASENAME, rasterized=True)

	# fix axis tick spacing
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(30))
	ax2.xaxis.set_major_locator(ticker.MultipleLocator(30))

	ax1.yaxis.set_major_locator(ticker.MultipleLocator(20))


# axis labels
ax1.set_ylabel(r'$\mathsf{Temperature}$' +r'$\mathsf{(K)}$', fontsize=12)
ax1.set_xlabel(r'$\mathsf{Latitude}$', fontsize=12)
ax2.set_ylabel(r'$\mathsf{Difference}$' + '\n' +r'$\mathsf{(K)}$', fontsize=12)
ax2.set_xlabel(r'$\mathsf{Latitude}$', fontsize=12)

plt.text(-0.10, 1.0, letters[n], fontsize=12, fontweight="bold", transform=ax1.transAxes)
n=n+1
plt.text(-0.10, 1.0, letters[n], fontsize=12, fontweight="bold", transform=ax2.transAxes)
n=n+1

#sort legend
handles, labels = ax1.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles2, labels2 = zip(*hl)


# Bottom section ---------------------------------------------------------------------

labs = [r'$\mathsf{Pressure}$' + '\n' +r'$\mathsf{(hPa)}$', '', r'$\mathsf{Potential}$' + '\n' +r'$\mathsf{Temperature}$'+r'$\mathsf{(K)}$', '', '', r'$\mathsf{Difference}$' + '\n' +r'$\mathsf{(K)}$',  r'$\mathsf{Pressure}$' + '\n' +r'$\mathsf{(hPa)}$', '', r'$\mathsf{Wet}$'+r'$\mathsf{Bulb}$'+ '\n'+r'$\mathsf{Potential}$' +'\n'+r'$\mathsf{Temperature}$'+r'$\mathsf{(K)}$', '', '', r'$\mathsf{Difference}$' + '\n' +r'$\mathsf{(K)}$', r'$\mathsf{Pressure}$' + '\n' +r'$\mathsf{(hPa)}$', '', r'$\mathsf{Humidity}$' + '\n' +r'$\mathsf{(g/kg)}$', '', '', r'$\mathsf{Difference}$' + '\n' +r'$\mathsf{(x10 g/kg)}$', r'$\mathsf{Pressure}$' + '\n' +r'$\mathsf{(hPa)}$', r'$\mathsf{Vertical Velocity}$'+'\n' + r'$\mathsf{(Pressure)}$' + '\n' +r'$\mathsf{(hPa/s)}$', '', r'$\mathsf{Difference}$' + '\n' +r'$\mathsf{(hPa/s)}$', '']

# keep track of which field/row we're on
row=0
# keep track of which gridspace/column we're plotting in for abs val
a = 0
# keep track of which gridspace/column we're plotting in for diff		
d = 0
# keep track of which vmin/max we're on
v = 0
# keep track of labels
l = 0
for p in fields:
	f = filenames[row]
	field = fields[row]

	# get out the data for the 1.0 S/So and 0.8 S/So
	presentcase = filebase + f + present +'.nc'
	eightcase = filebase + f + eight +'.nc'

	#plot the data -> S/S0 = 1.0
	ax = fig.add_subplot(absgrid[a])
	a=a+1

	if os.path.isfile(presentcase):
		ds = netCDF4.Dataset(presentcase)
		op = ds.variables[field][:]

		if n == 3:
			lat = ds.variables['latitude'][:]
		else:
			lat = ds.variables['lat'][:]
		p = ds.variables['lev'][:]
		ds.close() #close the file

		if n == 2:
			f = (1000/p)**0.286 # pressure
			t = np.squeeze(op)

			op = t*f[:,None]

		if n ==3:
			print('none')

		if n == 4:
			op = op*1000 # scale humidity up

		if n == 5:
			op = op*100 # scale vertical velocity 

		# color maps
		numcolor = op.size
		cm = plt.get_cmap(cmaps[row], 40)

		# create the colorbar
		if a > 9:
			tlevs=[-3.6, -3.2, -2.8, -2.4, -2.0,  -1.6,  -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2,  1.6,  2.0,  2.4,  2.8,  3.2, 3.6]
		else:
			tlevs = range(vmins[v], vmaxs[v])

		cnorm = colors.Normalize(vmins[v], vmaxs[v])
		scalarMap = mpl.cm.ScalarMappable(norm=cnorm, cmap=cm)
		
		# Create 2D lat/lon arrays for Basemap
		pressure2d, lat2d = np.meshgrid(p, lat)

		#ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in ts])
		ax.set_color_cycle([cm(1.*i/numcolor) for i in range(numcolor)])


		#plot the data
		c=plt.contourf(lat, p, np.squeeze(op), tlevs, cmap=cm, rasterized=True)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(30))

		# axis labels

		ax.set_ylabel(labs[l], fontsize=12)
		l = l+1

		ax.invert_yaxis()

		if n < 5:
			ax.set_xticklabels([])
			ax.set_xlabel(labs[l], fontsize=12)
			l = l+1

		if n == 5: ax.set_xlabel(r'$\mathsf{Latitude}$', fontsize=12)


		# This is the fix for the white lines between contour levels
		for i in c.collections:
  		  i.set_edgecolor("face")
		
		# add letter annotation
		plt.text(-0.15, 1.0, letters[n], fontsize=12, fontweight="bold", transform=ax.transAxes)

	#plot the data -> S/S0 = 0.8
	ax = fig.add_subplot(absgrid[a])
	a=a+1

	if os.path.isfile(eightcase):
		ds = netCDF4.Dataset(eightcase)
		o8 = ds.variables[field][:]
		if n == 3: 
			lat = ds.variables['latitude'][:]
		else:
			lat = ds.variables['lat'][:]
		p = ds.variables['lev'][:]
		ds.close() #close the file

		if n == 2:
			f = (1000/p)**0.286
			t = np.squeeze(o8)

			o8 = t*f[:,None]

		if n == 4:
			o8 = o8*1000

		if n == 5:
			o8 = o8*100

		#color based on temp
		numcolor = o8.size
		cm = plt.get_cmap(cmaps[row], 40)

		# create the colorbar
		if a > 9:
			tlevs=[-3.6, -3.2, -2.8, -2.4, -2.0,  -1.6,  -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2,  1.6,  2.0,  2.4,  2.8,  3.2, 3.6]
		else:
			tlevs = range(vmins[v], vmaxs[v])

		cnorm = colors.Normalize(vmins[v], vmaxs[v])
		scalarMap = mpl.cm.ScalarMappable(norm=cnorm, cmap=cm)
		
		# Create 2D lat/lon arrays for Basemap
		pressure2d, lat2d = np.meshgrid(p, lat)

		#ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in ts])
		ax.set_color_cycle([cm(1.*i/numcolor) for i in range(numcolor)])


		#plot the data
		c=plt.contourf(lat, p, np.squeeze(o8), tlevs, cmap=cm, rasterized=True)

		# axis labels
		ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
		ax.invert_yaxis()
		ax.set_yticklabels([])
	

		if n < 5:
			ax.set_xticklabels([])

		if n == 5: ax.set_xlabel(r'$\mathsf{Latitude}$', fontsize=12)


		# This is the fix for the white lines between contour levels
		for i in c.collections:
  		  i.set_edgecolor("face")

		# plot the colorbar - ABS value
		ax = fig.add_subplot(absgrid[a])
		a=a+1

		cb = plt.colorbar(c, cax=ax, label=labs[l])
		l = l+1

		tick_locator = ticker.MaxNLocator(nbins=3)
		cb.locator = tick_locator
		cb.update_ticks()
		v=v+1


	#plot the data - DIFF
	ax = fig.add_subplot(diffgrid[d])
	d=d+1
	if os.path.isfile(eightcase):
	
		#color based on temp
		numcolor = o8.size
		cm = plt.get_cmap(cmaps_d[row], 40)

		# create the colorbar
		if d == 7:
			tlevs=[-1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
		else:
			tlevs = range(vmins[v], vmaxs[v])

		cnorm = colors.Normalize(vmins[v], vmaxs[v])
		scalarMap = mpl.cm.ScalarMappable(norm=cnorm, cmap=cm)
		
		# Create 2D lat/lon arrays for Basemap
		pressure2d, lat2d = np.meshgrid(p, lat)

		#ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in ts])
		ax.set_color_cycle([cm(1.*i/numcolor) for i in range(numcolor)])

		if n == 4:
			diff =  (np.squeeze(o8) - np.squeeze(op))*10
		else:
			 diff = np.squeeze(o8) - np.squeeze(op)

		#plot the data
		c=plt.contourf(lat, p, diff, tlevs, cmap=cm, rasterized=True)

		# axis labels
		ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
		ax.invert_yaxis()
		ax.set_yticklabels([])


		ax.set_ylabel(labs[l], fontsize=12)
		l = l+1

		if n < 5:
			ax.set_xticklabels([])
			ax.set_xlabel(labs[l], fontsize=12)
			l = l+1

		if n == 5: ax.set_xlabel(r'$\mathsf{Latitude}$', fontsize=12)

		# This is the fix for the white lines between contour levels
		for i in c.collections:
  		  i.set_edgecolor("face")


		# plot the colorbar - DIFF value
		ax = fig.add_subplot(diffgrid[d])
		d=d+1
		cb = plt.colorbar(c, cax=ax)
		cb.set_label(label=labs[l], size=12)
		l = l+1

		tick_locator = ticker.MaxNLocator(nbins=3)
		cb.locator = tick_locator
		cb.update_ticks()
		v=v+1

	# go to next field/row
	n=n+1
	row = row+1


plt.show()

fig.savefig("figure3.pdf", bbox_inches='tight')

