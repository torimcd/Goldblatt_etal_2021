#!/usr/bin/env python3

import os
import sys
import numpy as np
import netCDF4
import operator
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D

sc = sys.argv[1]
run = 'Final'

# wet bulb potential temp
filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
outfilebase = 'wbpt_pressure_'

casename = sc

outfile = filebase + casename+'/'+outfilebase+'.nc'
print(outfile)

# check directly if the mergetime file exists
if os.path.isdir(filebase+casename):
	if not os.path.isfile(outfile):
		infile = filebase + casename+'/merged.nc'
		
		if not os.path.isfile(outfile):
			# calc potential temp
			syscall = 'cdo timmean -selyear,21/40 -select,name=T '+infile+' '+outfile
			print(syscall)
			os.system(syscall)



#create plot
#fig = plt.figure()
#ax = plt.subplot(211)

#dsloc = filebase + casename+'/'+outfilebase+'.nc'
#if os.path.isfile(dsloc):


	# open the merged file and get out the variables
#	ds = netCDF4.Dataset(dsloc)
#	t = ds.variables['T'][:]
#	lat = ds.variables['lat'][:]
#	p = ds.variables['lev'][:]
#	ds.close() #close the file


	#color based on temp
#	numcolor = t.size
#	cm = plt.get_cmap('Spectral_r', 20)

	# create the colorbar
#	tlevs = range(180, 325)

#	cnorm = colors.Normalize(180,325)
#	scalarMap = mpl.cm.ScalarMappable(norm=cnorm, cmap=cm)

	# Create 2D lat/lon arrays for Basemap
#	pressure2d, lat2d = np.meshgrid(p, lat)

	#ax.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in ts])
#	ax.set_color_cycle([cm(1.*i/numcolor) for i in range(numcolor)])


	#plot the data
#	c=plt.contourf(lat, p, np.squeeze(t), tlevs, cmap=cm)


	# This is the fix for the white lines between contour levels
#	for i in c.collections:
#  	  i.set_edgecolor("face")

#	cb = plt.colorbar(c, spacing='proportional')
#	cb.set_label('Temperature (K)')



#plt.title('Temperature - Solar Constant: '+str(round(int(sc)/1361,3)))
#plt.xlabel('Latitude')
#plt.ylabel('Pressure (hPa)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
#plt.gca().invert_yaxis()
#plt.tight_layout(w_pad=1.0, h_pad=5.0)

#plt.show()

#fig.savefig("zonaltemp_pressure_"+str(round(int(sc)/1361,3))+".pdf", bbox_inches='tight')

