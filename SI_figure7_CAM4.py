#!/usr/bin/env python3

import os
import sys
import numpy as np
from numpy.polynomial.polynomial import polyfit
import netCDF4
import operator
from matplotlib import colors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

run = 'Final'
years = 40

filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
outfilebase = 'anom_data'

casenames = {str(years)+'yr_1361': 'black', str(years)+'yr_1089': '#97dde5'}


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

lcl_c = []
z700_c = []
lts_c = []
qs850_c = []
tempsum_c = []
pal_c = []
eis_c = []
swcf_c = []
ah700_c = []
ah1000_c = []
ah_c = []
lcwp_c = []
cldlow_c = []
lhflx_c = []
ts_c = []

cldlow_anom = []
lhflx_anom = []
swcf_anom = []
ah_anom = []
eis_anom = []
lcwp_anom = []
ts_anom = []


for CASENAME in casenames.keys():
	if ((round(int(CASENAME[5:])/1361,3)) == 0.8):

		# GET LCL
		dsloc = filebase + CASENAME+'/'+outfilebase+'_lclv4.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcl_c = ds.variables['lcl'][:]
			ds.close() #close the file
			lcl_c = lcl_c.flatten()

		#GET z700
		dsloc = filebase + CASENAME+'/'+outfilebase+'_z700v4.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			z700_c = ds.variables['Z3'][:]
			ds.close() #close the file
			z700_c = z700_c.flatten()



		#GET lts
		dsloc = filebase + CASENAME+'/'+outfilebase+'ltsv4.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lts_c = ds.variables['lts'][:]
			ds.close() #close the file
			lts_c = lts_c.flatten()

		#GET qs850
		dsloc = filebase + CASENAME+'/'+outfilebase+'_qs850v3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			qs850_c = ds.variables['smr'][:]
			ds.close() #close the file
			qs850_c = qs850_c.flatten()	
	
		#GET tempsum
		dsloc = filebase + CASENAME+'/'+outfilebase+'_tempsumv3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			tempsum_c = ds.variables['T'][:]
			ds.close() #close the file
			tempsum_c = tempsum_c.flatten()

		#GET ah700
		dsloc = filebase + CASENAME+'/'+outfilebase+'_700ah3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah700_c = ds.variables['Q'][:]
			ds.close() #close the file
			ah700_c = ah700_c.flatten()

		#GET ah1000
		dsloc = filebase + CASENAME+'/'+outfilebase+'_1000ah3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah1000_c = ds.variables['QREFHT'][:]
			ds.close() #close the file
			ah1000_c = ah1000_c.flatten()


		#GET lcwp
		dsloc = filebase + CASENAME+'/'+outfilebase+'_lcwpv3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcwp_c = ds.variables['ICLDTWP'][:]
			ds.close() #close the file
			lcwp_c = lcwp_c.flatten()



		#GET swcf
		dsloc_swcf = filebase + CASENAME+'/'+outfilebase+'_swcfv3.nc'
		if os.path.isfile(dsloc_swcf):

			#print(dsloc_swcf)
			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_swcf)
			swcf_c = dscf.variables['SWCF'][:]
			dscf.close() #close the file
			swcf_c = swcf_c.flatten()


		#GET cldlow
		dsloc_cldlow = filebase + CASENAME+'/'+outfilebase+'_cldlowv3.nc'
		if os.path.isfile(dsloc_cldlow):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_cldlow)
			cldlow_c = dscf.variables['CLDLOW'][:]
			dscf.close() #close the file
			cldlow_c =cldlow_c.flatten()

		#GET lhflx
		dsloc_lhflx = filebase + CASENAME+'/'+outfilebase+'_lhflxv3.nc'
		if os.path.isfile(dsloc_lhflx):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_lhflx)
			lhflx_c = dscf.variables['LHFLX'][:]
			dscf.close() #close the file
			lhflx_c = lhflx_c.flatten()

		#GET ts
		dsloc = filebase + CASENAME+'/'+outfilebase+'_ts_v3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ts_c = ds.variables['TS'][:]
			ds.close() #close the file
			ts_c = ts_c.flatten()


		pal_c = (9.9)*(1-(1+2450000*qs850_c/(287.058*(tempsum_c/2)))/(1+(2450000**2)*qs850_c/(993*461.4*((tempsum_c/2)**2))))	

		eis_c = lts_c - pal_c*((z700_c/1000)-lcl_c)
		ah_c = list(map(operator.sub, ah1000_c*1000, ah700_c*1000))
		#mask = swcf_c < (-5)

	else:

		# GET LCL
		dsloc = filebase + CASENAME+'/'+outfilebase+'_lclv4.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcl = ds.variables['lcl'][:]
			ds.close() #close the file
			lcl = lcl.flatten()

		#GET z700
		dsloc = filebase + CASENAME+'/'+outfilebase+'_z700v4.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			z700 = ds.variables['Z3'][:]
			ds.close() #close the file
			z700 = z700.flatten()


		#GET lts
		dsloc = filebase + CASENAME+'/'+outfilebase+'ltsv4.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lts = ds.variables['lts'][:]
			ds.close() #close the file
			lts = lts.flatten()


		#GET qs850
		dsloc = filebase + CASENAME+'/'+outfilebase+'_qs850v3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			qs850 = ds.variables['smr'][:]
			ds.close() #close the file
			qs850 = qs850.flatten()	
	
		#GET tempsum
		dsloc = filebase + CASENAME+'/'+outfilebase+'_tempsumv3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			tempsum = ds.variables['T'][:]
			ds.close() #close the file
			tempsum = tempsum.flatten()

		#GET ah700
		dsloc = filebase + CASENAME+'/'+outfilebase+'_700ah3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah700 = ds.variables['Q'][:]
			ds.close() #close the file
			ah700 = ah700.flatten()

		#GET ah1000
		dsloc = filebase + CASENAME+'/'+outfilebase+'_1000ah3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			ah1000 = ds.variables['QREFHT'][:]
			ds.close() #close the file
			ah1000 = ah1000.flatten()


		#GET lcwp
		dsloc = filebase + CASENAME+'/'+outfilebase+'_lcwpv3.nc'
		if os.path.isfile(dsloc):

			# open the file and get out the variable
			ds = netCDF4.Dataset(dsloc)
			lcwp = ds.variables['ICLDTWP'][:]
			ds.close() #close the file
			lcwp = lcwp.flatten()



		#GET swcf
		dsloc_swcf = filebase + CASENAME+'/'+outfilebase+'_swcfv3.nc'
		if os.path.isfile(dsloc_swcf):

			#print(dsloc_swcf)
			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_swcf)
			swcf = dscf.variables['SWCF'][:]
			dscf.close() #close the file
			swcf = swcf.flatten()


		#GET cldlow
		dsloc_cldlow = filebase + CASENAME+'/'+outfilebase+'_cldlowv3.nc'
		if os.path.isfile(dsloc_cldlow):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_cldlow)
			cldlow = dscf.variables['CLDLOW'][:]
			dscf.close() #close the file
			cldlow =cldlow.flatten()

		#GET lhflx
		dsloc_lhflx = filebase + CASENAME+'/'+outfilebase+'_lhflxv3.nc'
		if os.path.isfile(dsloc_lhflx):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_lhflx)
			lhflx = dscf.variables['LHFLX'][:]
			dscf.close() #close the file
			lhflx = lhflx.flatten()


		#GET TS
		dsloc_ts = filebase + CASENAME+'/'+outfilebase+'_ts_v3.nc'
		if os.path.isfile(dsloc_ts):

			# open the file and get out the variable
			dscf = netCDF4.Dataset(dsloc_ts)
			ts = dscf.variables['TS'][:]
			dscf.close() #close the file
			ts = ts.flatten()


		pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	
		eis = lts - pal*((z700/1000)-lcl)
		#mask = swcf < (-5)
		ah = list(map(operator.sub, ah1000*1000, ah700*1000))


# eg eis_c is for sc=0.8, eis is sc=1.0 so eis_anom is 0.8sc - 1.0sc 
eis_anom = eis_c - eis
ah_anom = list(map(operator.sub, ah_c, ah))
lcwp_anom = lcwp_c - lcwp
cldlow_anom = cldlow_c - cldlow
swcf_anom = swcf_c -swcf
lhflx_anom = lhflx_c - lhflx
ts_anom = ts_c - ts

#fig = plt.figure()

#arr = plt.hist(ts_anom, bins=50)
#for i in range(50):
#    plt.text(arr[1][i],arr[0][i],str(arr[0][i]))

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

#ax1.scatter(eis_anom, swcf_anom, marker='o', alpha=0.2, color='orange')
#ax2.scatter(eis_anom, cldlow_anom, marker='o', alpha=0.2, color='skyblue')
#ax3.scatter(eis_anom, lcwp_anom, marker='o', alpha=0.2, color='yellowgreen')
#ax4.scatter(lhflx_anom, swcf_anom, marker='o', alpha=0.2, color='coral')
#ax5.scatter(lhflx_anom, cldlow_anom, marker='o', alpha=0.2, color='royalblue')
#ax6.scatter(lhflx_anom, lcwp_anom,marker='o', alpha=0.2, color='mediumseagreen')
#ax7.scatter(ah_anom, swcf_anom, marker='o', alpha=0.2, color='red')
#ax8.scatter(ah_anom, cldlow_anom, marker='o', alpha=0.2, color='blue')
#ax9.scatter(ah_anom, lcwp_anom, marker='o', alpha=0.2, color='darkgreen')

a = ax1.hist2d(lhflx_anom, cldlow_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.006, 0.002], [-0.00004, 0.00001]])
ax1.hlines(0, -0.1, 0.1)
ax1.vlines(0, -1, 1)

b = ax2.hist2d(ah_anom, cldlow_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0002], [-0.00003, 0.00001]])
ax2.hlines(0, -0.1, 0.1)
ax2.vlines(0, -1, 1)

c = ax3.hist2d(eis_anom, cldlow_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.005, 0.01], [-0.00003, 0.00001]])
ax3.hlines(0, -0.1, 0.1)
ax3.vlines(0, -1, 1)

d = ax4.hist2d(ts_anom, cldlow_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0004], [-0.00003, 0.00001]]) 
ax4.hlines(0, -0.1, 0.1)
ax4.vlines(0, -1, 1)

e = ax5.hist2d(lhflx_anom, lcwp_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.006, 0.002], [-0.02, 0.02]]) 
ax5.hlines(0, -0.1, 0.1)
ax5.vlines(0, -1, 1)

f = ax6.hist2d(ah_anom, lcwp_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0002], [-0.02, 0.02]])
ax6.hlines(0, -0.1, 0.1)
ax6.vlines(0, -1, 1)

g = ax7.hist2d(eis_anom, lcwp_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.005, 0.01], [-0.02, 0.02]])
ax7.hlines(0, -0.1, 0.1)
ax7.vlines(0, -1, 1)

h = ax8.hist2d(ts_anom, lcwp_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0004], [-0.02, 0.02]]) 
ax8.hlines(0, -0.1, 0.1)
ax8.vlines(0, -1, 1)

i = ax9.hist2d(lhflx_anom, swcf_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.006, 0.002], [-0.001, 0.006]])
ax9.hlines(0, -0.1, 0.1)
ax9.vlines(0, -1, 1)

j = ax10.hist2d(ah_anom, swcf_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0002], [-0.001, 0.006]])
ax10.hlines(0, -0.1, 0.1)
ax10.vlines(0, -1, 1)

k = ax11.hist2d(eis_anom, swcf_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.005, 0.01], [-0.001, 0.006]])
ax11.hlines(0, -0.1, 0.1)
ax11.vlines(0, -1, 1)

l = ax12.hist2d(ts_anom, swcf_anom, bins=50, cmap='Greys', norm=colors.Normalize(), vmin=0, vmax=100, range=[[-0.0004, 0.0004], [-0.001, 0.006]])
ax12.hlines(0, -0.1, 0.1)
ax12.vlines(0, -1, 1)

#ax1.set_title('EIS Anomaly vs Shortwave Cloud Forcing Anomaly')
#ax2.set_title('EIS Anomaly vs Low Cloud Fraction Anomaly')
#ax3.set_title('EIS Anomaly vs Low Cloud Water Path Anomaly')
#ax4.set_title('Latent Heat Flux Anomaly vs Shortwave Cloud Forcing Anomaly')
#ax5.set_title('Latent Heat Flux Anomaly vs Low Cloud Fraction Anomaly')
#ax6.set_title('Latent Heat Flux Anomaly vs Low Cloud Water Path Anomaly')
#ax7.set_title('Absolute Humidity Anomaly vs Shortwave Cloud Forcing Anomaly')
#ax8.set_title('Absolute Humidity vs Low Cloud Fraction Anomaly')
#ax9.set_title('Absolute Humidity vs Low Cloud Water Path Anomaly')

#plt.colorbar(a[3], ax=ax1)
#plt.colorbar(b[3], ax=ax2)
#plt.colorbar(c[3], ax=ax3)
#plt.colorbar(d[3], ax=ax4)
#plt.colorbar(e[3], ax=ax5)
#plt.colorbar(f[3], ax=ax6)
#plt.colorbar(g[3], ax=ax7)
#plt.colorbar(h[3], ax=ax8)
#plt.colorbar(i[3], ax=ax9)


ax1.set_ylabel('Low Cloud Fraction', fontsize=14)
#ax1.set_xlabel('Latent Heat Flux')
ax1.set_xticklabels([])

#ax2.set_ylabel('Low Cloud Fraction')
#ax2.set_xlabel('Absolute Humidity')

#ax3.set_ylabel('Low Cloud Water Path')
#ax3.set_xlabel('EIS')

#ax4.set_ylabel('SWCF')
#ax4.set_xlabel('Latent Heat Flux')

ax5.set_ylabel('Low Cloud Water Path', fontsize=14)
#ax5.set_xlabel('Latent Heat Flux')

#ax6.set_ylabel('Low Cloud Water Path')
#ax6.set_xlabel('Latent Heat Flux')

#ax7.set_ylabel('SWCF')
#ax7.set_xlabel('Humidity')

#ax8.set_ylabel('Low Cloud Fraction')
#ax8.set_xlabel('Humidity')

ax9.set_ylabel('SWCF', fontsize=14)
ax9.set_xlabel('Latent Heat Flux', fontsize=14)
ax9.invert_yaxis()

ax10.set_xlabel('Absolute Humidity', fontsize=14)
ax10.invert_yaxis()

ax11.set_xlabel('EIS', fontsize=14)
ax11.invert_yaxis()

ax12.set_xlabel('Surface Temperature', fontsize=14)
ax12.invert_yaxis()


# Best fit line
#plt.plot(np.unique(eis*round(int(CASENAME[5:])/1361,3)), np.poly1d(np.polyfit(eis*round(int(CASENAME[5:])/1361,3), mask*round(int(CASENAME[5:])/1361,3), 1))(np.unique(eis*round(int(CASENAME[5:])/1361,3))), linewidth=2)

# sort legend
#handles, labels = ax.get_legend_handles_labels()
#hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
#handles2, labels2 = zip(*hl)
# Display the legend below the plot
#plt.legend(handles2, labels2, bbox_to_anchor=(0., -0.75, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)

#plt.title('Shortwave Cloud Forcing vs Estimated Inversion Strength 60S to 60N')
#plt.xlabel('Estimated Inversion Strength (K)')
#plt.ylabel('Shortwave Cloud Forcing (W/m2)')
#ax.set_ylim([-120, 0])
#ax.set_xlim([-45, 45])


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

fig.savefig("anomaly_hist.pdf", bbox_inches='tight')


	

