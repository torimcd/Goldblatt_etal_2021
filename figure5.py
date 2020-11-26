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
from map_functions import get_map_field_anom
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap


# cloud forcing variables
lwcf = 'LWCF'
swcf = 'SWCF'
filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
#lwcfoutfilebase = 'lwcfyear_map_vl_diff_'
#swcfoutfilebase = 'swcfyear_map_vl_diff_'

outfilebase = 'map_annual_average'

casenames = ['1.1','1.05','1.0','0.95','0.9','0.85','0.8','0.75','0.7'] 

outer_grid = gridspec.GridSpec(2, 1, wspace=0.0, hspace=0.1)

#create plot
fig = plt.figure(figsize=(10, 8))
i=0
n=0
l=0

while i<2:
	inner_grid = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=outer_grid[i], wspace=0.0, hspace=0.15)	
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
		output = get_map_field_anom(field, CASENAME)

		fld = output[0]
		lats = output[1]
		lons = output[2]

	
		# setup the map
		m = Basemap(lat_0=0,lon_0=0,ax=ax)
		m.drawcoastlines()
		m.drawcountries()
		
		# Create 2D lat/lon arrays for Basemap
		lon2d, lat2d = np.meshgrid(lons, lats)
	
		# Plot 
		cs = m.pcolormesh(lon2d,lat2d,np.squeeze(fld), cmap='RdBu_r', latlon='True', vmin=-60, vmax=60)
		fig.add_subplot(ax)
	
		# This is the fix for the white lines between contour levels
		cs.set_edgecolor("face")
		
		# add letter annotation
		if l==0:
			plt.text(-0.10, 1.0, "a", fontsize=12, fontweight="bold", transform=ax.transAxes)
		elif l==1:
			plt.text(-0.10, 1.0, "b", fontsize=12, fontweight="bold", transform=ax.transAxes)

		l=l+2
	
		ax.set_title(CASENAME, fontsize=10)
	i=i+1

	
# Add Colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(cs, cax=cbar_ax, label=r'$\mathsf{W/m^2}$')

#plt.suptitle('Cloud Forcing Difference from Present')
plt.show()

fig.savefig("figure2.pdf", bbox_inches='tight')

