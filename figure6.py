#!/usr/bin/env python3
"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""
import matplotlib
matplotlib.use("Agg")

import os
import sys
import numpy as np
import netCDF4
import operator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap
import processing_functions as pf


# ------------------------------------------------------------------------
# change this section to match where you downloaded the model output files 
# ------------------------------------------------------------------------

download_path = '' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

filebase = download_path + 'FYSP_clouds_archive/CAM4/'
filebase_c5 = download_path + 'FYSP_clouds_archive/CAM5/'
outfileloc = download_path + 'temp_data/' # this is the location to save the processed netcdf files to

# ------------------------------------

casenames = ['07','0725','075','0775','08','0825','085','0875','09','0925','095','0975','10','1025','105','1075','11']
casenames_c5 = ['09','0925','095','0975','10','1025','105']

# field variables
fields = 'co2vmr,LHFLX,SHFLX,CLDLOW,CLDHGH,LWCF,SWCF,FSDS,FSDSC,FSNS,FSNSC,FLDS,FLDSC,FLNS,FLNSC'


pf.global_annual_average(filebase, outfileloc, fields, 'cam4')
pf.global_annual_average(filebase_c5, outfileloc, fields, 'cam5')
pf.prep_eis(filebase, outfileloc, 'cam4')
pf.prep_eis(filebase_c5, outfileloc, 'cam5')
pf.total_waterpath(filebase, outfileloc, 'cam4')


exp = [79219.8, 52021.4, 31909.4, 18299.8, 9552.6, 4780.9, 2154.7, 903.8, 368.9, 139.2, 51.4]
sc = [0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0,1.025,1.05]

sc_all = ['0.7','0.725','0.75','0.775','0.8','0.825','0.85','0.875','0.9','0.925','0.95','0.975','1.0','1.025','1.05','1.075','1.1']
sc_c5 =  ['0.9','0.925','0.95','0.975','1.0','1.025','1.05']

filenames = ['c4_global_average',
	'exp',
	'c4_global_average',
	'c4_global_average',
	'c4_eis',
	'c4_eis',
	'c4_global_average',
	'c4_global_average',
	'c4_wp_low',
	'c4_wp_high',
	'c4_global_average',
	'c4_global_average']

filenames_c5 = ['c5_global_average',
	'exp',
	'c5_global_average',
	'c5_global_average',
	'c5_eis',
	'c5_eis',
	'c5_global_average',
	'c5_global_average',
	'c5_wp_low',
	'c5_wp_high',
	'c5_global_average',
	'c5_global_average']

labels_c4 = ['CAM4',
	'expected',
	'Latent CAM4',
	'Sensible CAM4',
	'EIS CAM4',
	'eis',
	'Low CAM4',
	'High CAM4',
	'Low CAM4',
	'High CAM4',
	'Longwave CAM4',
	'Shortwave CAM4']

labels_c5 = ['CAM5',
	'expected',
	'Latent CAM5',
	'Sensible CAM5',
	'EIS CAM5',
	'eis',
	'Low CAM5',
	'High CAM5',
	'Low CAM5',
	'High CAM5',
	'Longwave CAM5',
	'Shortwave CAM5']

# field variables
fields= ['co2vmr',
	'exp',
	'LHFLX',
	'SHFLX',
	'EIS',
	'EIS',
	'CLDLOW',
	'CLDHGH',
	'ICLDTWP', 
	'ICLDTWP',
	'LWCF',
	'SWCF']

colors = ['black', 'black', 'saddlebrown', 'saddlebrown', 'peru', 'peru', 'blue', 'blue', 'blue', 'blue', 'lightskyblue', 'lightskyblue', 'blue', 'blue', 'lightskyblue', 'lightskyblue','salmon', 'salmon', 'red','red']

y_labels = [r'$\mathsf{CO_2}$'+ ' ' + r'$\mathsf{(ppm)}$',
	r'$\mathsf{Heat}$' + ' ' + r'$\mathsf{Flux}$' + ' ' +r'$\mathsf{(W/m^2)}$',
	'eis',
	r'$\mathsf{Cloud}$' + ' ' + r'$\mathsf{Fraction}$' ,
	r'$\mathsf{Cloud}$' + ' ' + r'$\mathsf{Water}$' + ' ' + r'$\mathsf{Path}$' + ' ' + r'$\mathsf{(W/m^2)}$',
	r'$\mathsf{Cloud}$' + ' ' + r'$\mathsf{Forcing}$' + ' ' +r'$\mathsf{(W/m^2)}$']

y_scales = [[0.1, 1000000],[0,100],[0,15],[0,0.5],[0,300],[-65,45]]

letters = ['a', 'b', 'c', 'd', 'e', 'f']


#create plot
fig = plt.figure(figsize=(7.08661,2.5))

# container 
grid = gridspec.GridSpec(2, 3, wspace=0.3, hspace=0.3)	

# keep track of which plot we're on
n=0
# keep track of which field is being plotted
f=0

# keep track of colors
c=0
for l in letters:
	ax = fig.add_subplot(grid[n])
	
	fn = filenames[f]
	fn_c5 = filenames_c5[f]
	field = fields[f]

	if (field == 'EIS'):
		i = 0
		eisplot = []
		eisplot_c5 = []

		ltsplot = []
		ltsplot_c5 = []
		
		for CASE in casenames:
			CASENAME = casenames[i]
			outfilebase = 'c4_eis_' + CASENAME + '_'
			lts_outfilebase = 'c4_lts_' + CASENAME + '_'	

			lcl = []
			z700 = []
			lts = []
			qs850 = []
			tempsum = []
			pal = []
			eis = []

	
			# GET LCL
			dsloc = outfileloc + outfilebase+'lcl.nc'
			if os.path.isfile(dsloc):

				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				lcl = ds.variables['lcl'][:]
				ds.close() #close the file
				lcl = lcl.flatten()
		
			#GET z700
			dsloc = outfileloc + outfilebase+'z700.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				z700 = ds.variables['Z3'][:]
				ds.close() #close the file
				z700 = z700.flatten()
		
		
			#GET lts
			dsloc = outfileloc + lts_outfilebase+'lts.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				lts = ds.variables['lts'][:]
				ds.close() #close the file
				lts = lts.flatten()
		
			#GET qs850
			dsloc = outfileloc + outfilebase+'qs850.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				qs850 = ds.variables['smr'][:]
				ds.close() #close the file
				qs850 = qs850.flatten()	
		
			#GET tempsum
			dsloc = outfileloc + outfilebase+'tempsum.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				tempsum = ds.variables['T'][:]
				ds.close() #close the file
				tempsum = tempsum.flatten()
		
			pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	
		
			eis = lts - pal*((z700/1000)-lcl)

			eisplot.append(eis.item(0))
			ltsplot.append(lts.item(0))
			i = i+1


		ax.plot(sc_all, eisplot, marker='o', ms=3, markeredgewidth=0.0, color='green', 
		label='EIS CAM4', rasterized=False)

		ax.plot(sc_all, ltsplot, marker='o', ms=3, markeredgewidth=0.0, color='lightgreen', 
		label='LTS CAM4', rasterized=False)


		
		i = 0
		for CASE in casenames_c5:
			CASENAME = casenames_c5[i]

			outfilebase = 'c5_eis_' + CASENAME + '_'
			lts_outfilebase = 'c5_lts_'	 + CASENAME + '_'	

			lcl = []
			z700 = []
			lts = []
			qs850 = []
			tempsum = []
			pal = []
			eisc5 = []

			# GET LCL
			dsloc = outfileloc + outfilebase +'lcl.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				lcl = ds.variables['lcl'][:]
				ds.close() #close the file
				lcl = lcl.flatten()
		
			#GET z700
			dsloc = outfileloc + outfilebase+'z700.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				z700 = ds.variables['Z3'][:]
				ds.close() #close the file
				z700 = z700.flatten()
		
		
			#GET lts
			dsloc = outfileloc + lts_outfilebase+'lts.nc'
			if os.path.isfile(dsloc):
		
				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				lts = ds.variables['lts'][:]
				ds.close() #close the file
				lts = lts.flatten()
		
			#GET qs850
			dsloc = outfileloc + outfilebase+'qs850.nc'
			if os.path.isfile(dsloc):

				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				qs850 = ds.variables['smr'][:]
				ds.close() #close the file
				qs850 = qs850.flatten()	
		
			#GET tempsum
			dsloc = outfileloc + outfilebase+'tempsum.nc'
			if os.path.isfile(dsloc):

				# open the file and get out the variable
				ds = netCDF4.Dataset(dsloc)
				tempsum = ds.variables['T'][:]
				ds.close() #close the file
				tempsum = tempsum.flatten()
		
			
			if os.path.isdir(outfileloc):
				pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	
		
				eisc5 = lts - pal*((z700/1000)-lcl)
				ltsplot_c5.append(lts.item(0))
				
				eisplot_c5.append(eisc5.item(0))
				i = i+1

		ax.plot(sc_c5, ltsplot_c5, marker='v', ms=3, markeredgewidth=0.0, linestyle='--', color='lightgreen', 
		label='LTS CAM5', rasterized=False)
		
		ax.plot(sc_c5, eisplot_c5, marker='v', ms=3, markeredgewidth=0.0, linestyle='--', color='green', 
		label='EIS CAM5', rasterized=False)



		ax.set_xlim([0.675,1.125])
		ax.set_ylabel(r'$\mathsf{Potential}$' + ' ' + r'$\mathsf{Temperature}$' + ' ' +r'$\mathsf{(K)}$', fontsize=5)
		ax.tick_params(labelsize=6) 

		plt.setp(ax.get_xticklabels(), fontsize=5)
		plt.setp(ax.get_yticklabels(), fontsize=5)
		

		


	else:
		# CAM4
		plotarray = []
		i=0
		for CASE in casenames:
			CASENAME = casenames[i]
			fp_c4 = outfileloc + fn + '_' + CASENAME +'.nc'
			
			if os.path.isfile(fp_c4):
	
				# open the file and get out the variable
				dsyear = netCDF4.Dataset(fp_c4)
				p = dsyear.variables[field][:]
				dsyear.close() #close the file
				p = p.flatten()

				if (field == 'co2vmr'):
					p=p*1000000

				plotarray.append(p.item(0))
			i +=1
	
		if len(plotarray) > 0:
			#plot the data
			ax.plot(sc_all, plotarray, marker='o', ms=3, markeredgewidth=0.0, color=colors[c], label=labels_c4[f], rasterized=False)
			ax.set_xlim([0.675,1.125])
			ax.set_ylabel(y_labels[n], fontsize=5)
			ax.tick_params(labelsize=6) 

		if n > 2:
			ax.set_xlabel(r'$\mathsf{S/S_0}$', fontsize=5)
		c +=1
	
		i=0
		plotarray_c5 = []
		for CASE in casenames_c5:
			CASENAME = casenames_c5[i]	
			fp_c5 = outfileloc + fn_c5 + '_' + CASENAME +'.nc'
			if os.path.isfile(fp_c5):
	
				# open the file and get out the variable
				dsyear = netCDF4.Dataset(fp_c5)
				p = dsyear.variables[field][:]
				dsyear.close() #close the file
				p = p.flatten()

				if (field == 'co2vmr'):
					p=p*1000000

				plotarray_c5.append(p.item(0))
				i +=1
	
		if len(plotarray_c5) > 0:
			#plot the data
			ax.plot(sc_c5, plotarray_c5, marker='v',  ms=3, markeredgewidth=0.0, linestyle='--', alpha=1, color=colors[c], label=labels_c5[f], rasterized=False)

		c +=1

	f = f+1
	fn = filenames[f]
	fn_c5 = filenames_c5[f]
	field = fields[f]


	if (field == 'exp'):
		ax.plot(sc, exp, marker='.', ms=3, color='lightgrey', label='Expectation', rasterized=False)

	else:
		# CAM4
		i=0
		plotarray = []
		for CASE in casenames:
			CASENAME = casenames[i]
			fp_c4 = outfileloc + fn + '_' + CASENAME + '.nc'
			if os.path.isfile(fp_c4):
	
				# open the file and get out the variable
				dsyear = netCDF4.Dataset(fp_c4)
				p = dsyear.variables[field][:]
				dsyear.close() #close the file
				p = p.flatten()

				if (field == 'co2vmr'):
					p=p*10000
	
				plotarray.append(p.item(0))
				i +=1


		if len(plotarray) > 0:
			#plot the data
			ax.plot(sc_all, plotarray, marker='o',  ms=3, markeredgewidth=0.0, color=colors[c], label=labels_c4[f], rasterized=False)

		c +=1
	
		i=0
		plotarray_c5 = []
		for CASE in casenames_c5:
			CASENAME = casenames_c5[i]	
			fp_c5 = outfileloc + fn_c5 + '_' + CASENAME +'.nc'
			if os.path.isfile(fp_c5):
	
				# open the file and get out the variable
				dsyear = netCDF4.Dataset(fp_c5)
				p = dsyear.variables[field][:]
				dsyear.close() #close the file
				p = p.flatten()

				if (field == 'co2vmr'):
					p=p*10000

				plotarray_c5.append(p.item(0))
				i +=1
	
		if len(plotarray_c5) > 0:
			#plot the data
			ax.plot(sc_c5, plotarray_c5, marker='v', ms=3, markeredgewidth=0.0, linestyle='--', alpha=1, color=colors[c], label=labels_c5[f], rasterized=False)

		c +=1

	# add letter annotation
	plt.text(-0.22, 1.0, letters[n], fontsize=6, fontweight="bold", transform=ax.transAxes)
	ax.set_ylim(y_scales[n])
	ax.tick_params(labelsize=5) 


	# Display the legend 
	plt.legend(loc=0, frameon=False, fontsize=5, prop={"size":5})
	if n == 0:
		ax.set_yscale("log")
	
	f = f+1
	n=n+1


plt.show()

fig.savefig("figures_main/figure6.pdf", format='pdf',bbox_inches='tight')

