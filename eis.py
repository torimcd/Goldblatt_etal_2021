#!/usr/bin/env python3

import os
import sys
import numpy
import netCDF4
import operator
import matplotlib.pyplot as plt

# This script calculates the Estimated Inversion Strength (EIS) out of various model output fields

filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
filebase_c5 = '/home/vmcd/projects/modelOutput/gcmexpt/CAM5/'

outfilebase = 'eis_'
lts_outfilebase = 'lts_'


casenames = {'0.7','0.725','0.75','0.775','0.8','0.825','0.85','0.875','0.9','0.925','0.95','0.975','1.0','1.025','1.05','1.075','1.1'}
casenames_c5 = {'0.9','0.925','0.95','0.975','1.0','1.025','1.05'}


# CAM4
for c in casenames:
	
	qs850 = filebase + c+'/'+outfilebase+'_qs850.nc'
	temp700 = filebase + c + '/'+outfilebase+'_temp700.nc'
	tempsurf = filebase + c+'/'+outfilebase+'_tempsurf.nc'
	tempsum = filebase + c+'/'+outfilebase+'_tempsum.nc'
	z700 = filebase + c + '/'+outfilebase+'_z700.nc'

	lcl = filebase + c+'/'+outfilebase+'_lcl.nc'

	# check directly if the mergetime file exists
	if os.path.isdir(filebase+c):
		if not os.path.isfile(lcl):
			infile = filebase + c+'/merged.nc'

			if not os.path.isfile(qs850):
				# calc saturation mixing ratio at 850 hPa
				syscall = 'cdo fldmean -timmean -selyear,21/40 -sellevidx,23 -select,name=smr -expr,\'smr=(Q/RELHUM)\'  '+infile+' '+qs850
				os.system(syscall)
		
			if not os.path.isfile(temp700):
				# get temp at 700hPa
				syscall = 'cdo fldmean -timmean -selyear,21/40 -sellevidx,21 -select,name=T '+infile+' '+temp700
				os.system(syscall)

			if not os.path.isfile(tempsurf):
				# get temp at surface
				syscall = 'cdo fldmean -timmean -selyear,21/40 -sellevidx,26 -select,name=T '+infile+' '+tempsurf
				os.system(syscall)

			if not os.path.isfile(z700):
				# get height of 700hPa
				syscall = 'cdo fldmean -timmean -selyear,21/40 -sellevidx,21 -select,name=Z3 '+infile+' '+z700
				os.system(syscall)

			if not os.path.isfile(lcl):
				# calc lifting condensation level
				syscall = 'cdo  fldmean -timmean -selyear,21/40 -sellevidx,26 -select,name=lcl -expr,\'lcl=((0.02+(TS-273.15)/5)*0.1*(1-RELHUM))\'  '+infile+' '+lcl
				os.system(syscall)

			if not os.path.isfile(tempsum):
				# temp sum 
				syscall = 'cdo add '+temp700+' '+tempsurf + ' '+tempsum
				os.system(syscall)


# calc temp CAM5
for c in casenames_c5:
	
	qs850_c5 = filebase_c5 + c+'/'+outfilebase+'_qs850.nc'
	temp700_c5 = filebase_c5 + c + '/'+outfilebase+'_temp700.nc'
	tempsurf_c5 = filebase_c5 + c+'/'+outfilebase+'_tempsurf.nc'
	tempsum_c5 = filebase_c5 + c+'/'+outfilebase+'_tempsum.nc'
	z700_c5 = filebase_c5 + c + '/'+outfilebase+'_z700.nc'

	lcl_c5 = filebase_c5 + c+'/'+outfilebase+'_lcl.nc'

	# check directly if the mergetime file exists
	if os.path.isdir(filebase_c5+c):
		if not os.path.isfile(lcl_c5):
			infile = filebase_c5 + c+'/merged.nc'

			if not os.path.isfile(qs850_c5):
				# calc saturation mixing ratio at 850 hPa
				syscall = 'cdo fldmean -timmean -selyear,31/60 -sellevidx,24 -select,name=smr -expr,\'smr=(Q/RELHUM)\'  '+infile+' '+qs850_c5
				os.system(syscall)
		
			if not os.path.isfile(temp700_c5):
				# get temp at 700hPa
				syscall = 'cdo fldmean -timmean -selyear,31/60 -sellevidx,21 -select,name=T '+infile+' '+temp700_c5
				os.system(syscall)

			if not os.path.isfile(tempsurf_c5):
				# get temp at surface
				syscall = 'cdo fldmean -timmean -selyear,31/60 -sellevidx,30 -select,name=T '+infile+' '+tempsurf_c5
				os.system(syscall)

			if not os.path.isfile(z700_c5):
				# get height of 700hPa
				syscall = 'cdo fldmean -timmean -selyear,31/60 -sellevidx,21 -select,name=Z3 '+infile+' '+z700_c5
				os.system(syscall)

			if not os.path.isfile(lcl_c5):
				# calc lifting condensation level
				syscall = 'cdo  fldmean -timmean -selyear,31/60 -sellevidx,30 -select,name=lcl -expr,\'lcl=((0.02+(TS-273.15)/5)*0.1*(1-RELHUM))\'  '+infile+' '+lcl_c5
				os.system(syscall)

			if not os.path.isfile(tempsum_c5):
				# temp sum 
				syscall = 'cdo add '+temp700_c5+' '+tempsurf_c5 + ' '+tempsum_c5
				os.system(syscall)

#create plot
fig = plt.figure()
ax = plt.subplot(211)

i = 0
eisplot_c4 = []
ltsplot_c4 = []
x_c4 = []
for c in casenames:

	lcl = []
	z700 = []
	lts = []
	qs850 = []
	tempsum = []
	pal = []
	eis = []


	
	# GET LCL
	dsloc = filebase + c+'/'+outfilebase+'_lcl.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		lcl = ds.variables['lcl'][:]
		ds.close() #close the file
		lcl = lcl.flatten()

	#GET z700
	dsloc = filebase + c+'/'+outfilebase+'_z700.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		z700 = ds.variables['Z3'][:]
		ds.close() #close the file
		z700 = z700.flatten()


	#GET lts
	dsloc = filebase + c+'/'+lts_outfilebase+'sub.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		lts = ds.variables['lts'][:]
		ds.close() #close the file
		lts = lts.flatten()
		ltsplot_c4.append(lts.item(0))

	#GET qs850
	dsloc = filebase + c+'/'+outfilebase+'_qs850.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		qs850 = ds.variables['smr'][:]
		ds.close() #close the file
		qs850 = qs850.flatten()	

	#GET tempsum
	dsloc = filebase + c+'/'+outfilebase+'_tempsum.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		tempsum = ds.variables['T'][:]
		ds.close() #close the file
		tempsum = tempsum.flatten()



	pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	

	eis = lts - pal*((z700/1000)-lcl)
	eisplot_c4.append(eis.item(0))
	x_c4.append(float(c))


	i = i+1


i = 0
eisplot_c5 = []
ltsplot_c5 = []
x_c5 = []
for c in casenames_c5:

	lcl = []
	z700 = []
	lts = []
	qs850 = []
	tempsum = []
	pal = []
	eis = []

	
	# GET LCL
	dsloc = filebase_c5 + c+'/'+outfilebase+'_lcl.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		lcl = ds.variables['lcl'][:]
		ds.close() #close the file
		lcl = lcl.flatten()

	#GET z700
	dsloc = filebase_c5 + c+'/'+outfilebase+'_z700.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		z700 = ds.variables['Z3'][:]
		ds.close() #close the file
		z700 = z700.flatten()


	#GET lts
	dsloc = filebase_c5 + c+'/'+lts_outfilebase+'sub.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		lts = ds.variables['lts'][:]
		ds.close() #close the file
		lts = lts.flatten()
		ltsplot_c5.append(lts.item(0))

	#GET qs850
	dsloc = filebase_c5 + c+'/'+outfilebase+'_qs850.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		qs850 = ds.variables['smr'][:]
		ds.close() #close the file
		qs850 = qs850.flatten()	

	#GET tempsum
	dsloc = filebase_c5 + c+'/'+outfilebase+'_tempsum.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		ds = netCDF4.Dataset(dsloc)
		tempsum = ds.variables['T'][:]
		ds.close() #close the file
		tempsum = tempsum.flatten()



	if os.path.isdir(filebase_c5+c):

		# equation 5 from Wood & Bretherton 2006
		pal = (9.9)*(1-(1+2450000*qs850/(287.058*(tempsum/2)))/(1+(2450000**2)*qs850/(993*461.4*((tempsum/2)**2))))	

		# equation 4 from Wood & Bretherton 2006
		eis = lts - pal*((z700/1000)-lcl)
		eisplot_c5.append(eis.item(0))
		x_c5.append(float(c))

	i = i+1

#sort the data
newx_c4, newy_eis = zip(*sorted(zip(x_c4,eisplot_c4)))
newx_c4, newy_lts = zip(*sorted(zip(x_c4,ltsplot_c4)))
newx_c5, c5_eis = zip(*sorted(zip(x_c5,eisplot_c5)))
newx_c5, c5_lts = zip(*sorted(zip(x_c5,ltsplot_c5)))


# make the plot
plt.plot(newx_c4, newy_lts, marker='o', linestyle='-', alpha = 0.7, color='black', label='LTS - CAM4')
plt.plot(newx_c5, c5_lts, marker='s', linestyle='--', alpha = 0.7, color='red', label='LTS - CAM5')
plt.plot(newx_c4, newy_eis, marker='o', linestyle='-', alpha = 0.7,color='b', label='EIS - CAM4')
plt.plot(newx_c5, c5_eis, marker='s', linestyle='--', alpha = 0.7, color='limegreen', label='EIS - CAM5')


plt.title('Global Average Lower Tropospheric Stability and Estimated Inversion Strength')
plt.xlabel('Solar Constant')
plt.ylabel('Potential Temperature (K)')
plt.axis([0.675,1.125,0,15])
plt.legend(loc=4)

plt.show()

fig.savefig("eis.pdf", bbox_inches='tight')

