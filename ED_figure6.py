#!/usr/bin/env python3
"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

This script prepares the model output to make a set of 2D histograms of cloud behavior against forcings. 
The actual figure is plotted in MATLAB, see /matlab_files/ED_figure6/

"""
import matplotlib
matplotlib.use("Agg")


import os
import sys
import numpy as np
import netCDF4
import operator
from matplotlib import colors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import processing_functions as pf

# ------------------------------------------------------------------------
# change this section to match where you downloaded the model output files 
# ------------------------------------------------------------------------

download_path = '' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

filebase = download_path + 'FYSP_clouds_archive/CAM4/'
filebase_c5 = download_path + 'FYSP_clouds_archive/CAM5/'

outfileloc = download_path + 'temp_data/' # this is the location to save the processed netcdf files to

pf.prep_anomaly_histograms(filebase, filebase_c5, outfileloc)

casenames_c4 = {'_10', '_09', '_08'}
casenames_c5 = {'_10', '_09'}


# Calc eis anom and plot v swcf
i = 0
lcl = []
z700 = []
lts = []
qs850 = []
tempsum = []
pal = []
eis = []
swcf = []
ah700 = []
ah1000 = []
ah = []
lcwp = []
cldlow = []
lhflx = []
ts = []

lcl_8 = []
z700_8 = []
lts_8 = []
qs850_8 = []
tempsum_8 = []
pal_8 = []
eis_8 = []
swcf_8 = []
ah700_8 = []
ah1000_8 = []
ah_8 = []
lcwp_8 = []
cldlow_8 = []
lhflx_8 = []
ts_8 = []

lcl_9 = []
z700_9 = []
lts_9 = []
qs850_9 = []
tempsum_9 = []
pal_9 = []
eis_9 = []
swcf_9 = []
ah700_9 = []
ah1000_9 = []
ah_9 = []
lcwp_9 = []
cldlow_9 = []
lhflx_9 = []
ts_9 = []

c4_cldlow_anom_8 = []
c4_lhflx_anom_8 = []
c4_swcf_anom_8 = []
c4_ah_anom_8 = []
c4_eis_anom_8 = []
c4_lcwp_anom_8 = []
c4_ts_anom_8 = []

c4_cldlow_anom_9 = []
c4_lhflx_anom_9 = []
c4_swcf_anom_9 = []
c4_ah_anom_9 = []
c4_eis_anom_9 = []
c4_lcwp_anom_9 = []
c4_ts_anom_9 = []

c5_cldlow_anom_9 = []
c5_lhflx_anom_9 = []
c5_swcf_anom_9 = []
c5_ah_anom_9 = []
c5_eis_anom_9 = []
c5_lcwp_anom_9 = []
c5_ts_anom_9 = []


for CASENAME in casenames_c4:
	if CASENAME == '_08':

		# GET LCL
		dsloc = outfileloc + 'c4_eis_map' + CASENAME + '_lcl.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcl_8 = ds.variables['lcl'][:]
			ds.close() #close the file
			lcl_8 = lcl_8.flatten()

		#GET z700
		dsloc = outfileloc + 'c4_eis_map' + CASENAME +'_z700.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			z700_8 = ds.variables['Z3'][:]
			ds.close() #close the file
			z700_8 = z700_8.flatten()


		#GET lts
		dsloc = outfileloc + 'c4_lts_map' + CASENAME +'_lts.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lts_8 = ds.variables['lts'][:]
			ds.close() #close the file
			lts_8 = lts_8.flatten()

		#GET qs850
		dsloc = outfileloc + 'c4_eis_map'+ CASENAME +'_qs850.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			qs850_8 = ds.variables['smr'][:]
			ds.close() #close the file
			qs850_8 = qs850_8.flatten()	
	
		#GET tempsum
		dsloc = outfileloc + 'c4_eis_map' + CASENAME +'_tempsum.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			tempsum_8 = ds.variables['T'][:]
			ds.close() #close the file
			tempsum_8 = tempsum_8.flatten()

		#GET ah700
		dsloc = outfileloc + 'c4_700ah' + CASENAME +'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah700_8 = ds.variables['Q'][:]
			ds.close() #close the file
			ah700_8 = ah700_8.flatten()

		#GET ah1000
		dsloc = outfileloc + 'c4_1000ah' + CASENAME +'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah1000_8 = ds.variables['QREFHT'][:]
			ds.close() #close the file
			ah1000_8 = ah1000_8.flatten()


		#GET lcwp
		dsloc = outfileloc + 'c4_map_wp_low' + CASENAME + '.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcwp_8 = ds.variables['ICLDTWP'][:]
			ds.close() #close the file
			lcwp_8 = lcwp_8.flatten()



		#GET swcf
		dsloc_swcf = outfileloc + 'c4_swcf'+CASENAME+'.nc'
		if os.path.isfile(dsloc_swcf):

			#print(dsloc_swcf)
			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_swcf)
			swcf_8 = dscf.variables['SWCF'][:]
			dscf.close() #close the file
			swcf_8 = swcf_8.flatten()


		#GET cldlow
		dsloc_cldlow = outfileloc + 'c4_cldlow' + CASENAME+'.nc'
		if os.path.isfile(dsloc_cldlow):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_cldlow)
			cldlow_8 = dscf.variables['CLDLOW'][:]
			dscf.close() #close the file
			cldlow_8 =cldlow_8.flatten()

		#GET lhflx
		dsloc_lhflx = outfileloc + 'c4_lhflx' + CASENAME+'.nc'
		if os.path.isfile(dsloc_lhflx):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_lhflx)
			lhflx_8 = dscf.variables['LHFLX'][:]
			dscf.close() #close the file
			lhflx_8 = lhflx_8.flatten()

		#GET ts
		dsloc = outfileloc + 'c4_ts' + CASENAME+'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ts_8 = ds.variables['TS'][:]
			ds.close() #close the file
			ts_8 = ts_8.flatten()


		pal_8 = (9.9)*(1-(1+2450000*qs850_8/(287.058*(tempsum_8/2)))/(1+(2450000**2)*qs850_8/(993*461.4*((tempsum_8/2)**2))))	

		eis_8 = lts_8 - pal_8*((z700_8/1000)-lcl_8)
		ah_8 = list(map(operator.sub, ah1000_8*1000, ah700_8*1000))


	elif CASENAME == '_09':
		# GET LCL
		dsloc = outfileloc + 'c4_eis_map' + CASENAME + '_lcl.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcl_9 = ds.variables['lcl'][:]
			ds.close() #close the file
			lcl_9 = lcl_9.flatten()

		#GET z700
		dsloc = outfileloc + 'c4_eis_map' + CASENAME +'_z700.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			z700_9 = ds.variables['Z3'][:]
			ds.close() #close the file
			z700_9 = z700_9.flatten()


		#GET lts
		dsloc = outfileloc + 'c4_lts_map' + CASENAME +'_lts.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lts_9 = ds.variables['lts'][:]
			ds.close() #close the file
			lts_9 = lts_9.flatten()

		#GET qs850
		dsloc = outfileloc + 'c4_eis_map'+ CASENAME +'_qs850.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			qs850_9 = ds.variables['smr'][:]
			ds.close() #close the file
			qs850_9 = qs850_9.flatten()	
	
		#GET tempsum
		dsloc = outfileloc + 'c4_eis_map' + CASENAME +'_tempsum.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			tempsum_9 = ds.variables['T'][:]
			ds.close() #close the file
			tempsum_9 = tempsum_9.flatten()

		#GET ah700
		dsloc = outfileloc + 'c4_700ah' + CASENAME +'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah700_9 = ds.variables['Q'][:]
			ds.close() #close the file
			ah700_9 = ah700_9.flatten()

		#GET ah1000
		dsloc = outfileloc + 'c4_1000ah' + CASENAME +'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah1000_9 = ds.variables['QREFHT'][:]
			ds.close() #close the file
			ah1000_9 = ah1000_9.flatten()


		#GET lcwp
		dsloc = outfileloc + 'c4_map_wp_low' + CASENAME + '.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcwp_9 = ds.variables['ICLDTWP'][:]
			ds.close() #close the file
			lcwp_9 = lcwp_9.flatten()



		#GET swcf
		dsloc_swcf = outfileloc + 'c4_swcf'+CASENAME+'.nc'
		if os.path.isfile(dsloc_swcf):

			#print(dsloc_swcf)
			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_swcf)
			swcf_9 = dscf.variables['SWCF'][:]
			dscf.close() #close the file
			swcf_9 = swcf_9.flatten()


		#GET cldlow
		dsloc_cldlow = outfileloc + 'c4_cldlow' + CASENAME+'.nc'
		if os.path.isfile(dsloc_cldlow):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_cldlow)
			cldlow_9 = dscf.variables['CLDLOW'][:]
			dscf.close() #close the file
			cldlow_9 =cldlow_9.flatten()

		#GET lhflx
		dsloc_lhflx = outfileloc + 'c4_lhflx' + CASENAME+'.nc'
		if os.path.isfile(dsloc_lhflx):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_lhflx)
			lhflx_9 = dscf.variables['LHFLX'][:]
			dscf.close() #close the file
			lhflx_9 = lhflx_9.flatten()

		#GET ts
		dsloc = outfileloc + 'c4_ts' + CASENAME+'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ts_9 = ds.variables['TS'][:]
			ds.close() #close the file
			ts_9 = ts_9.flatten()


		pal_9 = (9.9)*(1-(1+2450000*qs850_9/(287.058*(tempsum_9/2)))/(1+(2450000**2)*qs850_9/(993*461.4*((tempsum_9/2)**2))))	

		eis_9 = lts_9 - pal_9*((z700_9/1000)-lcl_9)
		ah_9 = list(map(operator.sub, ah1000_9*1000, ah700_9*1000))
		#mask = swcf_c < (-5)
		
	else:

		# GET LCL
		dsloc = outfileloc + 'c4_eis_map' + CASENAME + '_lcl.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcl = ds.variables['lcl'][:]
			ds.close() #close the file
			lcl = lcl.flatten()

		#GET z700
		dsloc = outfileloc + 'c4_eis_map' + CASENAME +'_z700.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			z700 = ds.variables['Z3'][:]
			ds.close() #close the file
			z700 = z700.flatten()


		#GET lts
		dsloc = outfileloc + 'c4_lts_map' + CASENAME +'_lts.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lts = ds.variables['lts'][:]
			ds.close() #close the file
			lts = lts.flatten()

		#GET qs850
		dsloc = outfileloc + 'c4_eis_map'+ CASENAME +'_qs850.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			qs850 = ds.variables['smr'][:]
			ds.close() #close the file
			qs850 = qs850.flatten()	
	
		#GET tempsum
		dsloc = outfileloc + 'c4_eis_map' + CASENAME +'_tempsum.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			tempsum = ds.variables['T'][:]
			ds.close() #close the file
			tempsum = tempsum.flatten()

		#GET ah700
		dsloc = outfileloc + 'c4_700ah' + CASENAME +'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah700 = ds.variables['Q'][:]
			ds.close() #close the file
			ah700 = ah700.flatten()

		#GET ah1000
		dsloc = outfileloc + 'c4_1000ah' + CASENAME +'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah1000 = ds.variables['QREFHT'][:]
			ds.close() #close the file
			ah1000 = ah1000.flatten()


		#GET lcwp
		dsloc = outfileloc + 'c4_map_wp_low' + CASENAME + '.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcwp = ds.variables['ICLDTWP'][:]
			ds.close() #close the file
			lcwp = lcwp.flatten()



		#GET swcf
		dsloc_swcf = outfileloc + 'c4_swcf'+CASENAME+'.nc'
		if os.path.isfile(dsloc_swcf):

			#print(dsloc_swcf)
			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_swcf)
			swcf = dscf.variables['SWCF'][:]
			dscf.close() #close the file
			swcf = swcf.flatten()


		#GET cldlow
		dsloc_cldlow = outfileloc + 'c4_cldlow' + CASENAME+'.nc'
		if os.path.isfile(dsloc_cldlow):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_cldlow)
			cldlow = dscf.variables['CLDLOW'][:]
			dscf.close() #close the file
			cldlow =cldlow.flatten()

		#GET lhflx
		dsloc_lhflx = outfileloc + 'c4_lhflx' + CASENAME+'.nc'
		if os.path.isfile(dsloc_lhflx):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_lhflx)
			lhflx = dscf.variables['LHFLX'][:]
			dscf.close() #close the file
			lhflx = lhflx.flatten()

		#GET ts
		dsloc = outfileloc + 'c4_ts' + CASENAME+'.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ts = ds.variables['TS'][:]
			ds.close() #close the file
			ts = ts.flatten()


		pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	

		eis = lts - pal*((z700/1000)-lcl)
		ah = list(map(operator.sub,ah1000*1000, ah700*1000))

# eis_anom is 0.8sc - 1.0sc 
c4_eis_anom_8 = list(map(operator.sub,eis_8, eis))
c4_ah_anom_8 = list(map(operator.sub, ah_8, ah))
c4_lcwp_anom_8 = list(map(operator.sub,lcwp_8, lcwp))
c4_cldlow_anom_8 = list(map(operator.sub, cldlow_8, cldlow))
c4_swcf_anom_8 = list(map(operator.sub,swcf_8, swcf))
c4_lhflx_anom_8 = list(map(operator.sub,lhflx_8, lhflx))
c4_ts_anom_8 = list(map(operator.sub,ts_8, ts))


#create plot
fig = plt.figure()

gs1 = gridspec.GridSpec(3,4)
gs1.update(wspace=0.05, hspace=0.05)

ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
ax3 = plt.subplot(gs1[2])
ax4 = plt.subplot(gs1[3])
ax5 = plt.subplot(gs1[4])
ax6 = plt.subplot(gs1[5])
ax7 = plt.subplot(gs1[6])
ax8 = plt.subplot(gs1[7])
ax9 = plt.subplot(gs1[8])
ax10 = plt.subplot(gs1[9])
ax11 = plt.subplot(gs1[10])
ax12 = plt.subplot(gs1[11])

a = ax1.hist2d(c4_lhflx_anom_8, c4_cldlow_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.006, 0.002], [-0.00004, 0.00001]])
ax1.hlines(0, -0.1, 0.1)
ax1.vlines(0, -1, 1)

b = ax2.hist2d(c4_ah_anom_8, c4_cldlow_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0002], [-0.00003, 0.00001]])
ax2.hlines(0, -0.1, 0.1)
ax2.vlines(0, -1, 1)

c = ax3.hist2d(c4_eis_anom_8, c4_cldlow_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.005, 0.01], [-0.00003, 0.00001]])
ax3.hlines(0, -0.1, 0.1)
ax3.vlines(0, -1, 1)

d = ax4.hist2d(c4_ts_anom_8, c4_cldlow_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0004], [-0.00003, 0.00001]]) 
ax4.hlines(0, -0.1, 0.1)
ax4.vlines(0, -1, 1)

e = ax5.hist2d(c4_lhflx_anom_8, c4_lcwp_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.006, 0.002], [-0.02, 0.02]]) 
ax5.hlines(0, -0.1, 0.1)
ax5.vlines(0, -1, 1)

f = ax6.hist2d(c4_ah_anom_8, c4_lcwp_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0002], [-0.02, 0.02]])
ax6.hlines(0, -0.1, 0.1)
ax6.vlines(0, -1, 1)

g = ax7.hist2d(c4_eis_anom_8, c4_lcwp_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.005, 0.01], [-0.02, 0.02]])
ax7.hlines(0, -0.1, 0.1)
ax7.vlines(0, -1, 1)

h = ax8.hist2d(c4_ts_anom_8, c4_lcwp_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0004], [-0.02, 0.02]]) 
ax8.hlines(0, -0.1, 0.1)
ax8.vlines(0, -1, 1)

i = ax9.hist2d(c4_lhflx_anom_8, c4_swcf_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.006, 0.002], [-0.001, 0.006]])
ax9.hlines(0, -0.1, 0.1)
ax9.vlines(0, -1, 1)

j = ax10.hist2d(c4_ah_anom_8, c4_swcf_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0002], [-0.001, 0.006]])
ax10.hlines(0, -0.1, 0.1)
ax10.vlines(0, -1, 1)

k = ax11.hist2d(c4_eis_anom_8, c4_swcf_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.005, 0.01], [-0.001, 0.006]])
ax11.hlines(0, -0.1, 0.1)
ax11.vlines(0, -1, 1)

l = ax12.hist2d(c4_ts_anom_8, c4_swcf_anom_8, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0004], [-0.001, 0.006]])
ax12.hlines(0, -0.1, 0.1)
ax12.vlines(0, -1, 1)



ax1.set_ylabel('Low Cloud Fraction', fontsize=6)
#ax1.set_xlabel('Latent Heat Flux')
ax1.set_xticklabels([])

ax5.set_ylabel('Low Cloud Water Path', fontsize=6)

ax9.set_ylabel('SWCF', fontsize=6)
ax9.set_xlabel('Latent Heat Flux', fontsize=6)
ax9.invert_yaxis()

ax10.set_xlabel('Absolute Humidity', fontsize=6)
ax10.invert_yaxis()

ax11.set_xlabel('EIS', fontsize=6)
ax11.invert_yaxis()

ax12.set_xlabel('Surface Temperature', fontsize=6)
ax12.invert_yaxis()


ax2.set_yticklabels([])
ax2.set_xticklabels([])

ax3.set_yticklabels([])
ax3.set_xticklabels([])

ax4.set_yticklabels([])
ax4.set_xticklabels([])

#ax.set_yticklabels([])
ax5.set_xticklabels([])

ax6.set_yticklabels([])
ax6.set_xticklabels([])

ax7.set_yticklabels([])
ax7.set_xticklabels([])

ax8.set_yticklabels([])
ax8.set_xticklabels([])

ax10.set_yticklabels([])

ax11.set_yticklabels([])
ax12.set_yticklabels([])

plt.setp(ax9.get_xticklabels(), rotation=45)
plt.setp(ax10.get_xticklabels(), rotation=45)
plt.setp(ax11.get_xticklabels(), rotation=45)
plt.setp(ax12.get_xticklabels(), rotation=45)

plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.015, 0.8])
plt.colorbar(a[3], cax = cax)

plt.show()

fig.savefig("figures_ED/ED_figure6.pdf", bbox_inches='tight')