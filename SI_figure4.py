#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

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

#run = 'Final'


# field variables
fields= ['T',
	'sigmaw',
	'Q',
	'OMEGA']

cmaps=['cool', 'cool','GnBu', 'RdGy']
cmaps_d=['RdBu_r','RdBu_r','BrBG','RdGy']

filenames = ['potentialtemp_pressure_',
	'wbpt_',
	'zonalhumidity_pressure_',
	'zonalvertvelc5_']

letters = ['a', 'b', 'c', 'd', 'e', 'f']

vmins = [200, -20, 240, -6, 0, -20,-5, -1.8] 
vmaxs = [602, 22, 320, 7, 20, 21, 6, 1.8]

filebase = '/nazko/home/shared/vmcd/modelOutput/gcmexpt/CAM5/'

present = '1.0'
eight = '0.9'

casenames = {'0.9': '#97dde5','0.925': '#2c7fb8','0.95': '#4c4cff','0.975':'#9366af','1.0': 'black','1.025': '#cb181d','1.05': '#fcae91',}

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
	ax1.plot(lat, ts_a, color=casenames[CASENAME], linewidth=1.5, label=CASENAME)
	ax2.plot(lat, ts_d, color=casenames[CASENAME], linewidth=1.5, label=CASENAME)

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

# Display the legend below the plot
#ax3.legend(handles2, labels2, bbox_to_anchor=(0, 0, 1, 1), loc=10, ncol=9, mode="expand", borderaxespad=0., frameon=False, fontsize=10)
#ax3.set_xticklabels([])
#ax3.set_yticklabels([])
#ax3.set_xticks([])
#ax3.set_yticks([])

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
	presentcase = filebase + present+'/'+f+'.nc'
	eightcase = filebase + eight +'/'+f+'.nc'

	if field == 'OMEGA':
		presentcase = '/home/vmcd/projects/output/zonalvertvel_1.0.nc'
		eightcase = '/home/vmcd/projects/output/zonalvertvel_0.9.nc'

	
	#plot the data - PRESENT
	ax = fig.add_subplot(absgrid[a])
	a=a+1

	print(presentcase)
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
			f = (1000/p)**0.286
			t = np.squeeze(op)

			op = t*f[:,None]

		if n ==3:
			print('none')

		if n == 4:
			op = op*1000

		if n == 5:
			op = op*100

		#color based on temp
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
		c=plt.contourf(lat, p, np.squeeze(op), tlevs, cmap=cm)
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

	#plot the data - EIGHT
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
		c=plt.contourf(lat, p, np.squeeze(o8), tlevs, cmap=cm)

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
		c=plt.contourf(lat, p, diff, tlevs, cmap=cm)

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

fig.savefig("figure3_cam5.pdf", bbox_inches='tight')
