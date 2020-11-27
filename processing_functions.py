#!/usr/bin/env python3
"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""

import os
import sys
import scipy
import numpy as np
import netCDF4

#-----------------------------------------------------------------------------------------#
# This file contains functions to process the model output for plotting in the figures
# accompanying Goldblatt et al 2020. Most of these functions use the Climate Data Operators
# (CDO) to time average, zonal average, perform mass weighting, etc. is used to prepare the 
# model output to be plotted on the maps in figure 2
#-----------------------------------------------------------------------------------------#

# this sets which cases to process. For speed to recreate the figures from the paper, by default this just
# uses the S/S0 = 1.0 and S/S0 = 0.8. If you would like to process results from all cases uncomment the 
# second set of lines but note that this will increase the processing time dramatically

cases_cam4 = {'08','10'}
cases_cam5 = {'09','10'}

#cases_cam4 = {'07','0725','075','0775','08','0825','085','0875','09','0925','095','0975','10','1025','105','1075','11'}
#cases_cam5 = {'09','0925','095','0975','10','1025','105'}

#-----------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on the maps in figure 2
#-----------------------------------------------------------------------------------------#
def map_annual_average(filebase, outloc, cam_version, fields):
	# This function calculates the annual average for all variables we're mapping in figs 2 and 4 that aren't on a pressure level

	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase = outloc + 'c4_map_annual_average'	
		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase = outloc + 'c5_map_annual_average'
		# the CAM5 cases
		casenames = cases_cam5
		# for CAM5 look at years 31-60
		selyear = '-selyear,31/60'	

	else:
		print('must supply a supported CAM version - cam4 or cam5')


	# average the fields for each case
	for c in casenames:
		outfile = outfilebase +'_'+ c +'.nc'

		# check directly if the input file exists
		if os.path.isdir(filebase):
			# only run if the output file does not exist
			if not os.path.isfile(outfile):
				infile = filebase + cam_version + '_' + c +'.nc'

				syscall = r"//usr//bin//cdo -timmean -yearmean " + selyear + " -select,name="+fields+" "+infile+" "+outfile
				print(syscall)
				os.system(syscall)



#-----------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on the maps in figure 2.
# It is essentially performing the same function as map_annual_average except that vertical
# velocity is on multiple pressure levels, so we need to choose one to plot it on the map.
# This code selects 700 hPa as the level to plot the vertical velocity
#-----------------------------------------------------------------------------------------#
def map_vert_velocity(filebase, outloc, cam_version):
	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase = outloc + 'c4_map_vert_velocity'	
		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase = outloc + 'c5_map_vert_velocity'
		# the CAM5 cases
		casenames = cases_cam5
		# for CAM5 look at years 31-60
		selyear = '-selyear,31/60'	

	else:
		print('must supply a supported CAM version - cam4 or cam5')

	for c in casenames:
		outfile = outfilebase +'_'+ c +'.nc'
	
		# check directly if the file exists
		if not os.path.isfile(outfile):
			if os.path.isdir(filebase):
				infile = filebase + cam_version + '_' + c +'.nc'
			
				syscall = r"//usr//bin//cdo timmean " + selyear + " -sellevidx,21 -select,name=OMEGA"+" "+infile+ " " +outfile
				print(syscall)
				os.system(syscall)



#-----------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on the maps in figure 4
# It is essentially performing the same function as map_annual_average except that total 
# water path is calculated for low and high clouds. This function outputs two netCDF files,
# one for the high waterpath, one for the low waterpath.
#-----------------------------------------------------------------------------------------#
def map_total_waterpath(filebase, outloc, cam_version):
	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filenames
		outfilebase_high = outloc + 'c4_map_wp_high'	
		outfilebase_low = outloc + 'c4_map_wp_low'	
		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	else:
		print('this function is only supported for cam4')

	for c in casenames:
		outfile_high = outfilebase_high +'_'+ c +'.nc'
		outfile_low = outfilebase_low +'_'+ c +'.nc'
	
		# check directly if the file exists
		if not os.path.isfile(outfile_high):
			if os.path.isdir(filebase):
				infile = filebase + cam_version + '_' + c +'.nc'
			
				syscall_high = r"//usr//bin//cdo -vertsum -timmean -sellevidx,1/17 " + selyear + " -select,name=ICLDTWP " +infile+ " " + outfile_high
				print(syscall_high)
				os.system(syscall_high)


		# check directly if the file exists
		if not os.path.isfile(outfile_low):
			if os.path.isdir(filebase):
				infile = filebase + cam_version + '_' + c +'.nc'

				syscall_low = r"//usr//bin//cdo -vertsum -timmean -sellevidx,21/26 " + selyear + " -select,name=ICLDTWP " +infile+ " " + outfile_low
				print(syscall_low)
				os.system(syscall_low)


#-----------------------------------------------------------------------------------------#
# This function is used to calculate Lower Tropospheric Stability (LTS) to be plotted on the 
# maps in figure 2. LTS is the difference in potential temperature between the suface and 
# 700 hPa, so we first select the temperature and pressure at 700 hPa, then calculate potential
# temperature at that level and at the surface, and finally subtract them.
#-----------------------------------------------------------------------------------------#
def prep_lts(filebase, outloc, cam_version):
	if cam_version == 'cam4':
		outfilebase = outloc + 'c4_lts_map'
		casenames = cases_cam4

		selyear = '-selyear,21/40'	
		plev = "696.79629"
	
	elif cam_version == 'cam5':
		outfilebase = outloc + 'c5_lts_map'
		casenames = cases_cam5

		selyear = '-selyear,31/60'	
		plev = "691.389430314302"

	else:
		print('must supply a supported CAM version - cam4 or cam5')


	for c in casenames:
		# set up the files to store the temporary output
		T700 = outfilebase+'_'+ c + '_700Temp.nc'
		outfile700 = outfilebase+'_'+ c + '_700.nc'
		outfile1000 = outfilebase+'_'+ c + '_1000.nc'
		outfilelts = outfilebase+'_'+ c + '_lts.nc'

		# check directly if the input file exists
		if os.path.isdir(filebase):
			if not os.path.isfile(outfilelts):
				infile = filebase + cam_version + '_' + c +'.nc'

				if not os.path.isfile(T700):
					# calc T at 700 hPa
					syscall = r"//usr//bin//cdo select,name=T -sellevel," + plev + " " +infile+ " " +T700
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(outfile700):
					# calc potential temp at 700hPa
					syscall = r"//usr//bin//cdo timmean " + selyear + " -select,name=lts -expr,\'lts=(T*(1000/" + plev + ")^0.286)\'  "+T700+" "+outfile700
					print(syscall)
					os.system(syscall)

			
				if not os.path.isfile(outfile1000):
					# calc potential temp at surface
					syscall = r"//usr//bin//cdo timmean " + selyear + " -select,name=lts -expr,\'lts=(TS*(1000/(PS*0.01))^0.286)\'  "+infile+" "+outfile1000
					print(syscall)
					os.system(syscall)

				# difference = lower tropospheric stability
				syscall = r"//usr//bin//cdo sub "+outfile700+" "+outfile1000 + " " + outfilelts
				print(syscall)
				os.system(syscall)



#-----------------------------------------------------------------------------------------#
# This function is used to calculate Estimated Inversion Strength (EIS) to be plotted on the 
# maps in figure 2. EIS is calculated using equation 4 from Wood & Bretherton 2006, which 
# requires calculating the variables necessary to calculate the moist potential temperature 
# gradient.
#-----------------------------------------------------------------------------------------#
def prep_eis(filebase, outloc, cam_version):
	# EIS uses LTS so run the LTS prep if not done already
	prep_lts(filebase, outloc, cam_version)

	if cam_version == 'cam4':
		outfilebase = outloc + 'c4_eis_map'
		casenames = cases_cam4

		selyear = '-selyear,21/40'	
	
	elif cam_version == 'cam5':
		outfilebase = outloc + 'c5_eis_map'
		casenames = cases_cam5

		selyear = '-selyear,31/60'	

	else:
		print('must supply a supported CAM version - cam4 or cam5')

	for c in casenames:
		qs850 = outfilebase+'_'+ c + '_qs850.nc'
		temp700 = outfilebase+'_'+ c + '_temp700.nc'
		tempsurf = outfilebase+'_'+ c + '_tempsurf.nc'
		tempsum = outfilebase+'_'+ c + '_tempsum.nc'
		z700 = outfilebase+'_'+ c + '_z700.nc'
		lcl = outfilebase+'_'+ c + '_lcl.nc'

		# check directly if the mergetime file exists
		if os.path.isdir(filebase):
			if not os.path.isfile(lcl):
				infile = filebase + cam_version + '_' + c +'.nc'

				if not os.path.isfile(qs850):
					# calc saturation mixing ratio at 850 hPa
					syscall = r"//usr//bin//cdo timmean " + selyear + " -sellevidx,23 -select,name=smr -expr,\'smr=(Q/RELHUM)\'  "+infile+" "+qs850
					print(syscall)
					os.system(syscall)
		
				if not os.path.isfile(temp700):
					# get temp at 700hPa
					syscall = r"//usr//bin//cdo timmean " + selyear + " -sellevidx,21 -select,name=T "+infile+" "+temp700
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(tempsurf):
					# get temp at surface
					syscall = r"//usr//bin//cdo timmean " + selyear + " -sellevidx,26 -select,name=T "+infile+" "+tempsurf
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(z700):
					# get height of 700hPa
					syscall = r"//usr//bin//cdo timmean " + selyear + " -sellevidx,21 -select,name=Z3 "+infile+" "+z700
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(lcl):
					# calc lifting condensation level
					syscall = r"//usr//bin//cdo  timmean " + selyear + " -sellevidx,26 -select,name=lcl -expr,\'lcl=((0.02+(TS-273.15)/5)*0.1*(1-RELHUM))\'  "+infile+" "+lcl
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(tempsum):
					# temp sum 
					syscall = r"//usr//bin//cdo add "+temp700+" "+tempsurf + " " +tempsum
					print(syscall)
					os.system(syscall)


#---------------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted zonally averaged in figure 3
#---------------------------------------------------------------------------------------------#
def zonal_average(filebase, outloc, am_version):
	# This function calculates the annual average for all variables we're mapping in figs 2 and 4 that aren't on a pressure level

	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase = outloc + 'c4_zonal_average'	
		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase = outloc + 'c5_zonal_average'
		# the CAM5 cases
		casenames = cases_cam5
		# for CAM5 look at years 31-60
		selyear = '-selyear,31/60'	

	else:
		print('must supply a supported CAM version - cam4 or cam5')


	# the fields we want to average for our plots  - these must not depend on pressure level
	fields = 'T, Q, OMEGA'

	# average the fields for each case
	for c in casenames:
		outfile = outfilebase +'_'+ c +'.nc'

		# check directly if the input file exists
		if os.path.isdir(filebase):
			# only run if the output file does not exist
			if not os.path.isfile(outfile):
				infile = filebase + cam_version + '_' + c +'.nc'

				syscall = r"//usr//bin//cdo zonmean -timmean "+ selyear + " -select,name="+ fields +" "+ infile+" "+outfile
				print(syscall)
				os.system(syscall)




#--------------------------------------------------------------------------#
# This function calculates the wet bulb potential temperature for figure 3
#-------------------------------------------------------------------------#
def wetbulb_potentialtemp(filebase, outloc, cam_version):

	# Wet bulb uses zonally averaged T so run the that first if not done already
	zonal_average(filebase, outloc, cam_version)

	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase =  outloc + 'c4_wetbulb'	
		zon_outfile =  outloc + 'c4_zonal_average'
		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase =  outloc + 'c5_wetbulb'
		zon_outfile =  outloc + 'c5_zonal_average'
		# the CAM5 cases
		casenames = cases_cam5
		# for CAM5 look at years 31-60
		selyear = '-selyear,31/60'	

	else:
		print('must supply a supported CAM version - cam4 or cam5')


	# calc wet bulb potential temp for each case
	for c in casenames:
		outfile = outfilebase +'_'+ c +'.nc' # the outfile for wet bulb
		zonfile = zon_outfile +'_'+ c +'.nc' # the zonal averaged T file to read

		# initialize these here so we can use them to write out the dimensions
		lat = []
		lev = []
		if os.path.isfile(zonfile):
			if not os.path.isfile(outfile):
				# open the zonal file and get out the variables
				ds = netCDF4.Dataset(dsloc_a)
				T = ds.variables['T'][:]
				lat = ds.variables['lat'][:]
				lev = ds.variables['lev'][:]
				ds.close() #close the file

				# set constants
				R= 287.1
				cp = 1004
				p0 = 1e5

				# unit conversion
				p = lev*100
				pp = np.ones(size(lat))*np.transpose(p)

				# calc potential temperature
				sigma = T*(p0/p)**(R/cp)

				# initialize
				sigmaw = 999*np.ones(size(sigma))

				# function we will integrate is
				# dTdp = pseudoadiabatig(p,T,Rd,Md,cpd,condensablegas,condensedphase)
				Rd = R
				cpd = cp
				Md = 0.02897
				rtol = 1**(-6)

				for i in sigmaw:
					pspan = np.array(pp[ii], 1e5)
					if pspan[0] == pspan[1]:
						sigmaw[i] = T[i]
					else:
						To = T[o]
						Tw = scipy.integrate.RK45(pseuedoadiabat(p, T, Rd, Md, cpd, 'h2o', 'l'), To, pspan, rtol=rtol)
						sigmaw[i] = Tw[-1]
		

		# save the calculates wbpt to a new netcdf file like the other processing steps
		try: ncfile.close()  # just to be safe, make sure dataset is not already open.
		except: pass
		ncfile = Dataset(outfile,mode='w',format='NETCDF4_CLASSIC') 

		lat_dim = ncfile.createDimension('lat', 73)     # latitude axis
		lev_dim = ncfile.createDimension('level', None) # unlimited axis (can be appended to).
		
		ncfile.title='WetBulb Potential Temperature'


		# Define two variables with the same names as dimensions,
		# a conventional way to define "coordinate variables".
		lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
		lat_var.units = 'degrees_north'
		lat_var.long_name = 'latitude'

		lev_var = ncfile.createVariable('lev', np.float64, ('lev',))
		lev_var.units = 'hPa'
		lev_var.long_name = 'pressure level'

		# Define a 2D variable to hold the data
		wbpt = ncfile.createVariable('sigmaw',np.float64,('lev','lat')) # note: unlimited dimension is leftmost
		wbpt.units = 'K' # degrees Kelvin
		wbpt.standard_name = 'wet_bulb_potential_temperature' # this is a CF standard name
		
		# write the values out
		lat_var[:] = lat
		lev_var[:] = lev
		wbpt[:,:] = sigmaw

		ncfile.close() # close the file when we're done writing it out


# function to integrate for calculating wet bulb potential temperature
def pseudoadiabatig():
	pass
    	

#-----------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on the maps in figure 5
#-----------------------------------------------------------------------------------------#
def cloud_forcing_all(filebase, outloc, fields, cam_version):

	cases_cam4 = {'07', '0725', '075', '0775','08', '0825', '085','0875','09', '0925','095', '0975', '10', '1025', '105','1075','11'}
	cases_cam5 = {'09', '0925','095', '0975', '10', '1025', '105'}
	
	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase =  outloc + 'c4_cloudforcing'	

		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase =  outloc + 'c5_cloudforcing'	

		# the CAM5 cases
		casenames = cases_cam5
		# for CAM5 look at years 31-60
		selyear = '-selyear,31/60'	

	else:
		print('must supply a supported CAM version - cam4 or cam5')
				

	# average the fields for each case
	for c in casenames:
		outfile = outfilebase +'_'+ c +'.nc'

		# check directly if the input file exists
		if os.path.isdir(filebase):
			# only run if the output file does not exist
			if not os.path.isfile(outfile):
				infile = filebase + cam_version + '_' + c +'.nc'

				syscall = r"//usr//bin//cdo -timmean -yearmean " + selyear + " -select,name="+fields+" "+infile+" "+outfile
				print(syscall)
				os.system(syscall)