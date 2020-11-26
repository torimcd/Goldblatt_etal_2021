#!/usr/bin/env python3

import os
import sys
import numpy
import netCDF4
import operator
import matplotlib.pyplot as plt

filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
filebase_c5 = '/home/vmcd/projects/modelOutput/gcmexpt/CAM5/'

outfilebase = 'lts_'


casenames = {'0.7','0.725','0.75','0.775','0.8','0.825','0.85','0.875','0.9','0.925','0.95','0.975','1.0','1.025','1.05','1.075','1.1'}
casenames_c5 = {'0.9','0.925','0.95','0.975','1.0','1.025','1.05'}



# CAM4
for c in casenames:
	T700 = filebase + c+'/'+outfilebase+'700Temp.nc'
	outfile700 = filebase + c+'/'+outfilebase+'700.nc'
	outfile1000 = filebase + c+'/'+outfilebase+'1000.nc'
	outfilesub = filebase + c+'/'+outfilebase+'sub.nc'

	# check directly if the mergetime file exists
	if os.path.isdir(filebase+c):
		if not os.path.isfile(outfilesub):
			infile = filebase + c+'/merged.nc'

			if not os.path.isfile(T700):
				# calc T at 700 hPa
				syscall = 'cdo sellevel,696.79629 -select,name=T '+infile+ ' '+T700
				os.system(syscall)

			if not os.path.isfile(outfile700):
				# calc potential temp 700
				syscall = 'cdo  fldmean -timmean -selyear,21/40 -select,name=lts -expr,\'lts=(T*(1000/696.79629)^0.286)\'  '+T700+' '+outfile700
				os.system(syscall)

			
			if not os.path.isfile(outfile1000):
				# calc potential temp at surface
				syscall = 'cdo fldmean -timmean -selyear,21/40 -select,name=lts -expr,\'lts=(TS*(1000/(PS*0.01))^0.286)\'  '+infile+' '+outfile1000
				os.system(syscall)

			# difference = lower tropospheric stability
			syscall = 'cdo sub '+outfile700+' '+outfile1000 + ' '+outfilesub
			os.system(syscall)


# calc temp CAM5
for c in casenames_c5:
	T700 = filebase_c5 + c+'/'+outfilebase+'_700Temp.nc'
	outfile700 = filebase_c5 + c+'/'+outfilebase+'700.nc'
	outfile1000 = filebase_c5 + c+'/'+outfilebase+'1000.nc'
	outfilesub = filebase_c5 + c+'/'+outfilebase+'sub.nc'

	# check directly if the mergetime file exists
	if os.path.isdir(filebase_c5+c):
		if not os.path.isfile(outfilesub):
			infile = filebase_c5 + c+'/merged.nc'

			if not os.path.isfile(T700):
				# calc T at 700 hPa
				syscall = 'cdo sellevel,691.389430314302 -select,name=T '+infile+ ' '+T700
				os.system(syscall)

			if not os.path.isfile(outfile700):
				# calc potential temp 700
				syscall = 'cdo  fldmean -timmean -selyear,31/60 -select,name=lts -expr,\'lts=(T*(1000/691.389430314302)^0.286)\'  '+T700+' '+outfile700
				os.system(syscall)

			
			if not os.path.isfile(outfile1000):
				# calc potential temp at surface
				syscall = 'cdo fldmean -timmean -selyear,31/60 -select,name=lts -expr,\'lts=(TS*(1000/(PS*0.01))^0.286)\'  '+infile+' '+outfile1000
				os.system(syscall)

			# difference = lower tropospheric stability
			syscall = 'cdo sub '+outfile700+' '+outfile1000 + ' '+outfilesub
			os.system(syscall)

#create plot
fig = plt.figure()
ax = plt.subplot(211)


x_c4=[]
y_c4=[]
for c in casenames:
	dsloc = filebase + c+'/'+outfilebase+'sub.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		dsyear = netCDF4.Dataset(dsloc)
		yr_lts = dsyear.variables['lts'][:]
		dsyear.close() #close the file
		yr_lts = yr_lts.flatten()

		# Save the new points for x and y
		x_c4.append(float(c))
		y_c4.append(yr_lts)




x_c5=[]
y_c5=[]
for c in casenames_c5:
	dsloc = filebase_c5 + c+'/'+outfilebase+'sub.nc'
	if os.path.isfile(dsloc):

		# open the file and get out the variable
		dsyear = netCDF4.Dataset(dsloc)
		yr_lts = dsyear.variables['lts'][:]
		dsyear.close() #close the file
		yr_lts = yr_lts.flatten()

		# Save the new points for x and y
		x_c5.append(float(c))
		y_c5.append(yr_lts)


#sort the data
newx, newy = zip(*sorted(zip(x_c4,y_c4)))
newx_c5, newy_c5 = zip(*sorted(zip(x_c5,y_c5)))

#plot it
plt.plot(newx, newy, marker='o', color='b', linestyle='--', label='CAM4')
plt.plot(newx_c5, newy_c5, marker='s', color='limegreen', linestyle='--', label='CAM5')



plt.title('Global Average Lower Tropospheric Stability')
plt.xlabel('Solar Constant')
plt.ylabel('Potential Temperature (K)')
plt.axis([0.675,1.125,7,14])
plt.legend(loc=4)

plt.show()

fig.savefig("lts.pdf", bbox_inches='tight')

