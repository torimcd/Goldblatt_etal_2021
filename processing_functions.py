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
import climt 
import json


#-----------------------------------------------------------------------------------------#
# This file contains functions to process the model output for plotting in the figures
# accompanying Goldblatt et al 2020. Most of these functions use the Climate Data Operators
# (CDO) to time average, zonal average, perform mass weighting, etc.
#-----------------------------------------------------------------------------------------#

# this sets which cases to process. For speed to recreate the figures from the paper, by default this just
# uses the S/S0 = 1.0 and S/S0 = 0.8. If you would like to process results from all cases uncomment the 
# second set of lines but note that this will increase the processing time dramatically

cases_cam4 = {'08','10'}
cases_cam5 = {'09','10'}

#cases_cam4 = {'07','0725','075','0775','08','0825','085','0875','09','0925','095','0975','10','1025','105','1075','11'}
#cases_cam5 = {'09','0925','095','0975','10','1025','105'}

#-----------------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on the maps in figures 2 and 5.
#-----------------------------------------------------------------------------------------------#
def map_annual_average(filebase, outloc, cam_version, fields):
	# This function calculates the annual average for all variables we're mapping in figs 2 and 5 that aren't on a pressure level

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
# This function is used to prepare the model output to be plotted on maps in figure 2.
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
# This function is used to prepare the model output to be plotted on maps in figure 5.
# It is essentially performing the same function as map_annual_average except that total 
# water path is calculated for low and high clouds. This function outputs two netCDF files,
# one for the high waterpath, one for the low waterpath.
#-----------------------------------------------------------------------------------------#
def map_total_waterpath(filebase, outloc, cam_version):
	cases_cam4 = {'08','09','10'} # we need the S/S0=0.9 case for these plots
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


#--------------------------------------------------------------------------------------------------------#
# This function is used to calculate Lower Tropospheric Stability (LTS) to be plotted on maps in figure 2.
# LTS is the difference in potential temperature between the suface and 700 hPa, so we first select the
# temperature and pressure at 700 hPa, then calculate potential temperature at that level and at the 
# surface, and finally subtract them.
#--------------------------------------------------------------------------------------------------------#
def map_prep_lts(filebase, outloc, cam_version):
	cases_cam4 = {'08','09','10'} # we need the S/S0=0.9 case for these plots

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



#--------------------------------------------------------------------------------------------------------#
# This function is used to calculate Estimated Inversion Strength (EIS) to be plotted on maps in figure 2. 
# EIS is calculated using equation 4 from Wood & Bretherton 2006, which requires calculating the 
# variables necessary to calculate the moist potential temperature gradient.
#--------------------------------------------------------------------------------------------------------#
def map_prep_eis(filebase, outloc, cam_version):
	cases_cam4 = {'08','09','10'} # we need the S/S0=0.9 case for these plots
	
	# EIS uses LTS so run the LTS prep if not done already
	map_prep_lts(filebase, outloc, cam_version)

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
# This function is used to prepare the model output to be plotted zonally averaged in figure 3.
#---------------------------------------------------------------------------------------------#
def zonal_average(filebase, outloc, cam_version, numfields):
	# This function calculates the annual average for all variables we're mapping that aren't on a pressure level
	if numfields == 'all':
		cases_cam4 = {'07','0725','075','0775', '08','0825','085','0875','09','0925','095','0975','10', '1025', '105', '1075','11'}
		cases_cam5 = {'09','0925','095','0975','10', '1025', '105'}

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
	fields = 'TS,T,Q,OMEGA'

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
# This function calculates the wet bulb potential temperature for figure 3.
#-------------------------------------------------------------------------#
def wetbulb_potentialtemp(filebase, outloc, cam_version):

	# prepare the variables to use in wet bulb calculation
	wbpt_prep(filebase, outloc, cam_version)

#--------------------------------------------------------------------------#
# This function prepares the variables needed to calculate wet bulb potential
# temperature for figure 3.
#-------------------------------------------------------------------------#
def wbpt_prep(filebase, outloc, cam_version):
	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase = outloc + 'c4_wbpt_prep'	
		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase = outloc + 'c5_wbpt_prep'
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

				syscall = r"//usr//bin//cdo -timmean "+ selyear + " -select,name=T "+ infile+" "+outfile
				print(syscall)
				os.system(syscall)




#-----------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on maps in figure 6.
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



#------------------------------------------------------------------------------------#
# This function is used to prepare the model output for figure 4.
# Calculates a global annual average of the provided fields for all cases.
#------------------------------------------------------------------------------------#
def global_annual_average(filebase, outloc, fields, cam_version):
	cases_cam4 = {'07', '0725', '075', '0775','08', '0825', '085','0875','09', '0925','095', '0975', '10', '1025', '105','1075','11'}
	cases_cam5 = {'09', '0925','095', '0975', '10', '1025', '105'}
	
	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filename 
		outfilebase =  outloc + 'c4_global_average'	

		# the CAM4 cases
		casenames = cases_cam4
		# for CAM4 look at years 21-40
		selyear = '-selyear,21/40'	

	elif cam_version == 'cam5':
		# the processed filename 
		outfilebase =  outloc + 'c5_global_average'	

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

				syscall = r"//usr//bin//cdo -timmean -fldmean -yearmean " + selyear + " -select,name="+fields+" "+infile+" "+outfile
				print(syscall)
				os.system(syscall)



#-----------------------------------------------------------------------------------------#
# This function is used to calculate Lower Tropospheric Stability (LTS) for figure 4. 
# LTS is the difference in potential temperature between the suface and 700 hPa, 
# so we first select the temperature and pressure at 700 hPa, then calculate potential
# temperature at that level and at the surface, and finally subtract them.
#-----------------------------------------------------------------------------------------#
def prep_lts(filebase, outloc, cam_version):
	cases_cam4 = {'07', '0725', '075', '0775','08', '0825', '085','0875','09', '0925','095', '0975', '10', '1025', '105','1075','11'}
	cases_cam5 = {'09', '0925','095', '0975', '10', '1025', '105'}
	
	if cam_version == 'cam4':
		outfilebase = outloc + 'c4_lts'
		casenames = cases_cam4

		selyear = '-selyear,21/40'	
		plev = "696.79629"
	
	elif cam_version == 'cam5':
		outfilebase = outloc + 'c5_lts'
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
					syscall = r"//usr//bin//cdo fldmean -timmean " + selyear + " -select,name=lts -expr,\'lts=(T*(1000/" + plev + ")^0.286)\'  "+T700+" "+outfile700
					print(syscall)
					os.system(syscall)

			
				if not os.path.isfile(outfile1000):
					# calc potential temp at surface
					syscall = r"//usr//bin//cdo fldmean -timmean " + selyear + " -select,name=lts -expr,\'lts=(TS*(1000/(PS*0.01))^0.286)\'  "+infile+" "+outfile1000
					print(syscall)
					os.system(syscall)

				# difference = lower tropospheric stability
				syscall = r"//usr//bin//cdo sub "+outfile700+" "+outfile1000 + " " + outfilelts
				print(syscall)
				os.system(syscall)



#-----------------------------------------------------------------------------------------#
# This function is used to calculate Estimated Inversion Strength (EIS) for figure 4. 
# EIS is calculated using equation 4 from Wood & Bretherton 2006, which requires calculating
# the variables necessary to calculate the moist potential temperature gradient.
#-----------------------------------------------------------------------------------------#
def prep_eis(filebase, outloc, cam_version):
	# EIS uses LTS so run the LTS prep if not done already
	prep_lts(filebase, outloc, cam_version)

	cases_cam4 = {'07', '0725', '075', '0775','08', '0825', '085','0875','09', '0925','095', '0975', '10', '1025', '105','1075','11'}
	cases_cam5 = {'09', '0925','095', '0975', '10', '1025', '105'}

	if cam_version == 'cam4':
		outfilebase = outloc + 'c4_eis'
		casenames = cases_cam4

		selyear = '-selyear,21/40'	
	
	elif cam_version == 'cam5':
		outfilebase = outloc + 'c5_eis'
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
					syscall = r"//usr//bin//cdo fldmean -timmean " + selyear + " -sellevidx,23 -select,name=smr -expr,\'smr=(Q/RELHUM)\'  "+infile+" "+qs850
					print(syscall)
					os.system(syscall)
		
				if not os.path.isfile(temp700):
					# get temp at 700hPa
					syscall = r"//usr//bin//cdo fldmean -timmean " + selyear + " -sellevidx,21 -select,name=T "+infile+" "+temp700
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(tempsurf):
					# get temp at surface
					syscall = r"//usr//bin//cdo fldmean -timmean " + selyear + " -sellevidx,26 -select,name=T "+infile+" "+tempsurf
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(z700):
					# get height of 700hPa
					syscall = r"//usr//bin//cdo fldmean -timmean " + selyear + " -sellevidx,21 -select,name=Z3 "+infile+" "+z700
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(lcl):
					# calc lifting condensation level
					syscall = r"//usr//bin//cdo  fldmean -timmean " + selyear + " -sellevidx,26 -select,name=lcl -expr,\'lcl=((0.02+(TS-273.15)/5)*0.1*(1-RELHUM))\'  "+infile+" "+lcl
					print(syscall)
					os.system(syscall)

				if not os.path.isfile(tempsum):
					# temp sum 
					syscall = r"//usr//bin//cdo add "+temp700+" "+tempsurf + " " +tempsum
					print(syscall)
					os.system(syscall)


#-----------------------------------------------------------------------------------------#
# This function is used to prepare the model output to be plotted on the maps in figure 5
# It is essentially performing the same function as map_annual_average except that total 
# water path is calculated for low and high clouds. This function outputs two netCDF files,
# one for the high waterpath, one for the low waterpath.
#-----------------------------------------------------------------------------------------#
def total_waterpath(filebase, outloc, cam_version):
	cases_cam4 = {'07', '0725', '075', '0775','08', '0825', '085','0875','09', '0925','095', '0975', '10', '1025', '105','1075','11'}
	cases_cam5 = {'09', '0925','095', '0975', '10', '1025', '105'}

	# check whether we're processessing CAM4 or CAM5
	if cam_version == 'cam4':
		# the processed filenames
		outfilebase_high = outloc + 'c4_wp_high'	
		outfilebase_low = outloc + 'c4_wp_low'	
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
			
				syscall_high = r"//usr//bin//cdo -vertsum -fldmean -timmean -sellevidx,1/17 " + selyear + " -select,name=ICLDTWP " +infile+ " " + outfile_high
				print(syscall_high)
				os.system(syscall_high)


		# check directly if the file exists
		if not os.path.isfile(outfile_low):
			if os.path.isdir(filebase):
				infile = filebase + cam_version + '_' + c +'.nc'

				syscall_low = r"//usr//bin//cdo -vertsum -fldmean -timmean -sellevidx,21/26 " + selyear + " -select,name=ICLDTWP " +infile+ " " + outfile_low
				print(syscall_low)
				os.system(syscall_low)


#-----------------------------------------------------------------------------------------#
# This function is used to generate the radiative fluxes from CAM5 (fluxes_*.txt files
# which are provided in data archive /FYSP_clouds_archive/radiative_transfer/CAM5/), 
# for Extended Data Figure 8.
#-----------------------------------------------------------------------------------------#
def generate_fluxes_cam5(filebase):
	
	# read in pressure/temp profile
	infile = filebase +'rad_tr_profile.txt'

	dtype1 = np.dtype([('Pressure','f8'),('Temperature','f8'),('H2O','f8'),('O3','f8'),('P_mid','f8'),('T_mid','f8'),('Water_mid','f8'),('Ozone_mid','f8')])
	profile = np.loadtxt(infile, dtype=dtype1, skiprows=1)

	p_int = np.array(profile['Pressure'])
	T_int = np.array(profile['Temperature'])
	Water = np.array(profile['H2O'])
	Ozone = np.array(profile['O3'])
	p_mid = np.array(profile['P_mid'])
	T_mid = np.array(profile['T_mid'])
	Water_mid = np.array(profile['Water_mid'])
	Ozone_mid = np.array(profile['Ozone_mid'])

	# Remove dummy data from last row of mid levels - dummy data added because np.loadtxt requires columns to be same length
	T_mid = np.resize(T_mid, 41)
	p_mid = np.resize(p_mid, 41)
	Water_mid = np.resize(Water_mid, 41)
	Ozone_mid = np.resize(Ozone_mid, 41)


	
	swradiation = climt.RRTMGShortwave()
	lwradiation = climt.RRTMGLongwave()
	surface = climt.SlabSurface()

	# Create shortwave state
	m_levs=dict(label='levels',values=np.arange(41), units='')
	i_levs=dict(label='levels', values=np.arange(42), units='')
	state=climt.get_default_state([swradiation, lwradiation, surface], mid_levels=m_levs, interface_levels=i_levs)

	# Change values to match my atmosphere profile
	state['surface_temperature'] = state['surface_temperature'] - 11.76
	state['air_pressure'].values = p_mid.reshape(1,1,41)
	state['air_pressure_on_interface_levels'].values = p_int.reshape(1,1,42)
	state['air_temperature'].values = T_mid.reshape(1,1,41)
	state['specific_humidity'].values = Water_mid.reshape(1,1,41)
	state['mole_fraction_of_ozone_in_air'].values = Ozone_mid.reshape(1,1,41)
	state['surface_albedo_for_direct_shortwave'] = state['surface_albedo_for_direct_shortwave'] +0.06
	state['single_scattering_albedo_due_to_cloud'] = state['single_scattering_albedo_due_to_cloud'] +0.06
	state['surface_albedo_for_diffuse_shortwave'] = state['surface_albedo_for_diffuse_shortwave'] +0.06
	state['single_scattering_albedo_due_to_aerosol'] = state['single_scattering_albedo_due_to_aerosol'] +0.06
	state['surface_albedo_for_diffuse_near_infrared'] = state['surface_albedo_for_diffuse_near_infrared'] +0.06
	state['surface_albedo_for_direct_near_infrared'] = state['surface_albedo_for_direct_near_infrared'] +0.06
	state['zenith_angle'] = state['zenith_angle'] + np.pi/3

	state['mole_fraction_of_methane_in_air'] = state['mole_fraction_of_methane_in_air'] + 0.00000184278
	state['mole_fraction_of_nitrous_oxide_in_air'] = state['mole_fraction_of_nitrous_oxide_in_air'] + 3.28*10**-7

	state['air_pressure'].attrs['units'] = 'hPa'
	state['air_pressure_on_interface_levels'].attrs['units'] = 'hPa'


	# See the state of all variables
	#print(state)

	# Calc unperturbed fluxes - CO2 at 330 ppm default
	swtendencies, swdiagnostics = swradiation(state)
	lwtendencies, lwdiagnostics = lwradiation(state)

	unpert_swup = np.array(swdiagnostics['upwelling_shortwave_flux_in_air']) 
	unpert_swdown = np.array(swdiagnostics['downwelling_shortwave_flux_in_air']) 
	unpert_lwup = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
	unpert_lwdown = np.array(lwdiagnostics['downwelling_longwave_flux_in_air']) 


	# prepare data to save to file
	swup_330_data = np.reshape(unpert_swup,42)*0.965529262*0.501472063
	swdown_330_data = np.reshape(unpert_swdown, 42)*0.965529262*0.501472063
	lwup_330_data = np.reshape(unpert_lwup,42)
	lwdown_330_data = np.reshape(unpert_lwdown,42)
	
	net_flux_330 = swdown_330_data + lwdown_330_data - swup_330_data - lwup_330_data
	net_flux_330_data = np.reshape(net_flux_330,42)


	# prepare data to save to file
	c5_pressure = np.reshape(state['air_pressure_on_interface_levels'].values, 42)

	if not os.path.isfile(filebase + "fluxes_330.txt"):
		# write the data out as json
		with open(filebase + "fluxes_330.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_330_data, 'sw_dn':swdown_330_data, 'lw_up':lwup_330_data, 'lw_dn':lwdown_330_data}, default=lambda x: list(x), indent=4))
	

		
	if not os.path.isfile(filebase + "fluxes_100.txt"):
		# ----------------------  Change CO2 concentration - 100ppm --------------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*0.3030303

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_100 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_100 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_100 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_100 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

		# prepare data to save to file
		swup_100_data = np.reshape(swup_100,42)*0.965529262*0.501472063
		swdown_100_data = np.reshape(swdown_100, 42)*0.965529262*0.501472063
		lwup_100_data = np.reshape(lwup_100,42)
		lwdown_100_data = np.reshape(lwdown_100,42)

		net_flux_100 = swdown_100_data + lwdown_100_data - swup_100_data - lwup_100_data -  net_flux_330
		net_flux_100_data = np.reshape(net_flux_100,42)

		c5_pressure = np.reshape(state['air_pressure_on_interface_levels'].values, 42)


		# write the data out as json
		with open(filebase + "fluxes_100.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_100_data, 'sw_dn':swdown_100_data, 'lw_up':lwup_100_data, 'lw_dn':lwdown_100_data}, default=lambda x: list(x), indent=4))



	if not os.path.isfile(filebase + "fluxes_500.txt"):
		# ------------------ Change CO2 concentration - 500ppm ------------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*1.51515152

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_500 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_500 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_500 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_500 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_500_data = np.reshape(swup_500,42)*0.965529262*0.501472063
		swdown_500_data = np.reshape(swdown_500, 42)*0.965529262*0.501472063
		lwup_500_data = np.reshape(lwup_500,42)
		lwdown_500_data = np.reshape(lwdown_500,42)

		net_flux_500 = swdown_500_data + lwdown_500_data - swup_500_data - lwup_500_data - net_flux_330
		net_flux_500_data = np.reshape(net_flux_500,42)

		# write it out as json
		with open(filebase + "fluxes_500.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_500_data, 'sw_dn':swdown_500_data, 'lw_up':lwup_500_data, 'lw_dn':lwdown_500_data}, default=lambda x: list(x), indent=4))



	if not os.path.isfile(filebase + "fluxes_1000.txt"):
		# --------------------------- Change CO2 concentration - 1000ppm ----------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*3.03030303

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_1000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_1000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_1000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_1000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_1000_data = np.reshape(swup_1000,42)*0.965529262*0.501472063
		swdown_1000_data = np.reshape(swdown_1000, 42)*0.965529262*0.501472063
		lwup_1000_data = np.reshape(lwup_1000,42)
		lwdown_1000_data = np.reshape(lwdown_1000,42)

		net_flux_1000 = swdown_1000_data + lwdown_1000_data - swup_1000_data - lwup_1000_data - net_flux_330
		net_flux_1000_data = np.reshape(net_flux_1000,42)

		# write the data out as json
		with open(filebase + "fluxes_1000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_1000_data, 'sw_dn':swdown_1000_data, 'lw_up':lwup_1000_data, 'lw_dn':lwdown_1000_data}, default=lambda x: list(x), indent=4))




	if not os.path.isfile(filebase + "fluxes_1500.txt"):
		# --------------------- Change CO2 concentration - 1500ppm -------------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*4.54545455

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_1500 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_1500 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_1500 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_1500 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_1500_data = np.reshape(swup_1500,42)*0.965529262*0.501472063
		swdown_1500_data = np.reshape(swdown_1500, 42)*0.965529262*0.501472063
		lwup_1500_data = np.reshape(lwup_1500,42)
		lwdown_1500_data = np.reshape(lwdown_1500,42)

		net_flux_1500 = swdown_1500_data + lwdown_1500_data - swup_1500_data - lwup_1500_data - net_flux_330
		net_flux_1500_data = np.reshape(net_flux_1500,42)

		with open(filebase + "fluxes_1500.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_1500_data, 'sw_dn':swdown_1500_data, 'lw_up':lwup_1500_data, 'lw_dn':lwdown_1500_data}, default=lambda x: list(x), indent=4))





	if not os.path.isfile(filebase + "fluxes_2000.txt"):
		# ------------------- Change CO2 concentration - 2000ppm -----------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*6.06060606

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_2000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_2000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_2000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_2000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_2000_data = np.reshape(swup_2000,42)*0.965529262*0.501472063
		swdown_2000_data = np.reshape(swdown_2000, 42)*0.965529262*0.501472063
		lwup_2000_data = np.reshape(lwup_2000,42)
		lwdown_2000_data = np.reshape(lwdown_2000,42)

		net_flux_2000 = swdown_2000_data + lwdown_2000_data - swup_2000_data - lwup_2000_data - net_flux_330
		net_flux_2000_data = np.reshape(net_flux_2000,42)

		with open(filebase + "fluxes_2000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_2000_data, 'sw_dn':swdown_2000_data, 'lw_up':lwup_2000_data, 'lw_dn':lwdown_2000_data}, default=lambda x: list(x), indent=4))




	if not os.path.isfile(filebase + "fluxes_5000.txt"):
		# ------------------------- Change CO2 concentration - 5000ppm -----------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*15.1515152

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_5000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_5000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_5000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_5000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_5000_data = np.reshape(swup_5000,42)*0.965529262*0.501472063
		swdown_5000_data = np.reshape(swdown_5000, 42)*0.965529262*0.501472063
		lwup_5000_data = np.reshape(lwup_5000,42)
		lwdown_5000_data = np.reshape(lwdown_5000,42)

		net_flux_5000 = swdown_5000_data + lwdown_5000_data - swup_5000_data - lwup_5000_data - net_flux_330
		net_flux_5000_data = np.reshape(net_flux_5000,42)


		with open(filebase + "fluxes_5000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_5000_data, 'sw_dn':swdown_5000_data, 'lw_up':lwup_5000_data, 'lw_dn':lwdown_5000_data}, default=lambda x: list(x), indent=4))




	if not os.path.isfile(filebase + "fluxes_10000.txt"):
		#  -------------------- Change CO2 concentration - 10000ppm --------------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*30.3030303

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_10000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_10000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_10000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_10000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_10000_data = np.reshape(swup_10000,42)*0.965529262*0.501472063
		swdown_10000_data = np.reshape(swdown_10000, 42)*0.965529262*0.501472063
		lwup_10000_data = np.reshape(lwup_10000,42)
		lwdown_10000_data = np.reshape(lwdown_10000,42)

		net_flux_10000 = swdown_10000_data + lwdown_10000_data - swup_10000_data - lwup_10000_data - net_flux_330
		net_flux_10000_data = np.reshape(net_flux_10000,42)


		with open(filebase + "fluxes_10000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_10000_data, 'sw_dn':swdown_10000_data, 'lw_up':lwup_10000_data, 'lw_dn':lwdown_10000_data}, default=lambda x: list(x), indent=4))



	if not os.path.isfile(filebase + "fluxes_20000.txt"):
		# ------------------- Change CO2 concentration - 20000ppm --------------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*60.6060606

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_20000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_20000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_20000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_20000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])



		# prepare data to save to file
		swup_20000_data = np.reshape(swup_20000,42)*0.965529262*0.501472063
		swdown_20000_data = np.reshape(swdown_20000, 42)*0.965529262*0.501472063
		lwup_20000_data = np.reshape(lwup_20000,42)
		lwdown_20000_data = np.reshape(lwdown_20000,42)

		net_flux_20000 = swdown_20000_data + lwdown_20000_data - swup_20000_data - lwup_20000_data - net_flux_330
		net_flux_20000_data = np.reshape(net_flux_20000,42)

		with open(filebase + "fluxes_20000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_20000_data, 'sw_dn':swdown_20000_data, 'lw_up':lwup_20000_data, 'lw_dn':lwdown_20000_data}, default=lambda x: list(x), indent=4))



	if not os.path.isfile(filebase + "fluxes_30000.txt"):
		# -------------------- Change CO2 concentration - 30000ppm -------------------------------------
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*90.9090909

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_30000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_30000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_30000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_30000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])



		# prepare data to save to file
		swup_30000_data = np.reshape(swup_30000,42)*0.965529262*0.501472063
		swdown_30000_data = np.reshape(swdown_30000, 42)*0.965529262*0.501472063
		lwup_30000_data = np.reshape(lwup_30000,42)
		lwdown_30000_data = np.reshape(lwdown_30000,42)

		net_flux_30000 = swdown_30000_data + lwdown_30000_data - swup_30000_data - lwup_30000_data - net_flux_330
		net_flux_30000_data = np.reshape(net_flux_30000,42)

		with open(filebase + "fluxes_30000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_30000_data, 'sw_dn':swdown_30000_data, 'lw_up':lwup_30000_data, 'lw_dn':lwdown_30000_data}, default=lambda x: list(x), indent=4))





	if not os.path.isfile(filebase + "fluxes_35000.txt"):
		# Change CO2 concentration - 35000ppm
		state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*106.060606

		# Calc perturbed fluxes
		swtendencies, swdiagnostics = swradiation(state)
		lwtendencies, lwdiagnostics = lwradiation(state)

		swup_35000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
		swdown_35000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])
		lwup_35000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
		lwdown_35000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])


		# prepare data to save to file
		swup_35000_data = np.reshape(swup_35000,42)*0.965529262*0.501472063
		swdown_35000_data = np.reshape(swdown_35000, 42)*0.965529262*0.501472063
		lwup_35000_data = np.reshape(lwup_35000,42)
		lwdown_35000_data = np.reshape(lwdown_35000,42)

		net_flux_35000 = swdown_35000_data + lwdown_35000_data - swup_35000_data - lwup_35000_data - net_flux_330
		net_flux_35000_data = np.reshape(net_flux_35000,42)

		with open(filebase + "fluxes_35000.txt", 'w+') as datafile_id:
			datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_35000_data, 'sw_dn':swdown_35000_data, 'lw_up':lwup_35000_data, 'lw_dn':lwdown_35000_data}, default=lambda x: list(x), indent=4))
	

# -------------------------------------------------------------------------------
# This function calculates the net surface, 200hPa, and top of atmosphere (TOA)
# radiative forcing for CAM5. CAM4, and the SMART line by line model, for Extended
# Data Figure 8.
# -------------------------------------------------------------------------------
def calculate_rad_fluxes(filebase):
	cam5_filebase = filebase + 'CAM5/'
	cam3_filebase = filebase + 'CAM3/'
	smart_filebase = filebase + 'SMART/'

	generate_fluxes_cam5(cam5_filebase)

	#--------------------------------------------------------
	# CAM5
	#--------------------------------------------------------

	# Read in fluxes
	# ------ 330
	cam5_infile = cam5_filebase + 'fluxes_330.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_330 = np.asarray(ds['lw_up'])
	c5lwdn_330 = np.asarray(ds['lw_dn'])
	c5swup_330 = np.asarray(ds['sw_up'])
	c5swdn_330 = np.asarray(ds['sw_dn'])

	c5_total_flux_330 = c5swup_330*(-1) + c5swdn_330 - c5lwup_330 + c5lwdn_330


	# ------ 100
	cam5_infile = cam5_filebase + 'fluxes_100.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_100 = np.asarray(ds['lw_up'])
	c5lwdn_100 = np.asarray(ds['lw_dn'])
	c5swup_100 = np.asarray(ds['sw_up'])
	c5swdn_100 = np.asarray(ds['sw_dn'])

	c5_total_flux_100 = c5swup_100*(-1) + c5swdn_100 - c5lwup_100 + c5lwdn_100 - c5_total_flux_330



	# ------ 500
	cam5_infile = cam5_filebase + 'fluxes_500.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_500 = np.asarray(ds['lw_up'])
	c5lwdn_500 = np.asarray(ds['lw_dn'])
	c5swup_500 = np.asarray(ds['sw_up'])
	c5swdn_500 = np.asarray(ds['sw_dn'])

	c5_total_flux_500 = c5swup_500*(-1) + c5swdn_500 - c5lwup_500 + c5lwdn_500 - c5_total_flux_330


	# ------ 1000
	cam5_infile = cam5_filebase + 'fluxes_1000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_1000 = np.asarray(ds['lw_up'])
	c5lwdn_1000 = np.asarray(ds['lw_dn'])
	c5swup_1000 = np.asarray(ds['sw_up'])
	c5swdn_1000 = np.asarray(ds['sw_dn'])

	c5_total_flux_1000 = c5swup_1000*(-1) + c5swdn_1000 - c5lwup_1000 + c5lwdn_1000 - c5_total_flux_330


	# ------ 1500
	cam5_infile = cam5_filebase + 'fluxes_1500.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_1500 = np.asarray(ds['lw_up'])
	c5lwdn_1500 = np.asarray(ds['lw_dn'])
	c5swup_1500 = np.asarray(ds['sw_up'])
	c5swdn_1500 = np.asarray(ds['sw_dn'])

	c5_total_flux_1500 = c5swup_1500*(-1) + c5swdn_1500 - c5lwup_1500 + c5lwdn_1500 - c5_total_flux_330


	# ------ 2000
	cam5_infile = cam5_filebase + 'fluxes_2000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_2000 = np.asarray(ds['lw_up'])
	c5lwdn_2000 = np.asarray(ds['lw_dn'])
	c5swup_2000 = np.asarray(ds['sw_up'])
	c5swdn_2000 = np.asarray(ds['sw_dn'])

	c5_total_flux_2000 = c5swup_2000*(-1) + c5swdn_2000 - c5lwup_2000 + c5lwdn_2000 - c5_total_flux_330


	# ------ 5000
	cam5_infile = cam5_filebase + 'fluxes_5000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_5000 = np.asarray(ds['lw_up'])
	c5lwdn_5000 = np.asarray(ds['lw_dn'])
	c5swup_5000 = np.asarray(ds['sw_up'])
	c5swdn_5000 = np.asarray(ds['sw_dn'])

	c5_total_flux_5000 = c5swup_5000*(-1) + c5swdn_5000 - c5lwup_5000 + c5lwdn_5000 - c5_total_flux_330


	# ------ 10000
	cam5_infile = cam5_filebase + 'fluxes_10000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_10000 = np.asarray(ds['lw_up'])
	c5lwdn_10000 = np.asarray(ds['lw_dn'])
	c5swup_10000 = np.asarray(ds['sw_up'])
	c5swdn_10000 = np.asarray(ds['sw_dn'])

	c5_total_flux_10000 = c5swup_10000*(-1) + c5swdn_10000 - c5lwup_10000 + c5lwdn_10000 - c5_total_flux_330


	# ------ 20000
	cam5_infile = cam5_filebase + 'fluxes_20000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_20000 = np.asarray(ds['lw_up'])
	c5lwdn_20000 = np.asarray(ds['lw_dn'])
	c5swup_20000 = np.asarray(ds['sw_up'])
	c5swdn_20000 = np.asarray(ds['sw_dn'])

	c5_total_flux_20000 = c5swup_20000*(-1) + c5swdn_20000 - c5lwup_20000 + c5lwdn_20000 - c5_total_flux_330


	# ------ 30000
	cam5_infile = cam5_filebase + 'fluxes_30000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_30000 = np.asarray(ds['lw_up'])
	c5lwdn_30000 = np.asarray(ds['lw_dn'])
	c5swup_30000 = np.asarray(ds['sw_up'])
	c5swdn_30000 = np.asarray(ds['sw_dn'])

	c5_total_flux_30000 = c5swup_30000*(-1) + c5swdn_30000 - c5lwup_30000 + c5lwdn_30000 - c5_total_flux_330


	# ------ 35000
	cam5_infile = cam5_filebase + 'fluxes_35000.txt'
	ds = json.load(open(cam5_infile))

	c5lwup_35000 = np.asarray(ds['lw_up'])
	c5lwdn_35000 = np.asarray(ds['lw_dn'])
	c5swup_35000 = np.asarray(ds['sw_up'])
	c5swdn_35000 = np.asarray(ds['sw_dn'])

	c5_total_flux_35000 = c5swup_35000*(-1) + c5swdn_35000 - c5lwup_35000 + c5lwdn_35000 - c5_total_flux_330





	#--------------------------------------------------------
	# CAM3
	#--------------------------------------------------------

	# ------ 330
	cam3_infile_330 = cam3_filebase + 'CliMT_CAM3_1D_330.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_330)
	
	cam3_flux_330_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:]
	cam3_flux_330_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19]
	cam3_flux_330_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:]

	cam3_pressure = ds.variables['Pres'][:]

	# ------ 100
	cam3_infile_100 = cam3_filebase + 'CliMT_CAM3_1D_100.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_100)
	
	cam3_flux_100_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_100_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_100_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa

	# ------ 500
	cam3_infile_500 = cam3_filebase + 'CliMT_CAM3_1D_500.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_500)
	
	cam3_flux_500_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_500_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_500_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa

	# ------ 1000
	cam3_infile_1000 = cam3_filebase + 'CliMT_CAM3_1D_1000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_1000)
	
	cam3_flux_1000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_1000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_1000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa


	# ------ 2000
	cam3_infile_2000 = cam3_filebase + 'CliMT_CAM3_1D_2000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_2000)
	
	cam3_flux_2000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_2000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_2000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa


	# ------ 5000
	cam3_infile_5000 = cam3_filebase + 'CliMT_CAM3_1D_5000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_5000)
	
	cam3_flux_5000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_5000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_5000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa


	# ------ 10000
	cam3_infile_10000 = cam3_filebase + 'CliMT_CAM3_1D_10000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_10000)
	
	cam3_flux_10000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_10000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_10000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa


	# ------ 20000
	cam3_infile_20000 = cam3_filebase + 'CliMT_CAM3_1D_20000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_20000)
	
	cam3_flux_20000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_20000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_20000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa


	# ------ 30000
	cam3_infile_30000 = cam3_filebase + 'CliMT_CAM3_1D_30000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_30000)
	
	cam3_flux_30000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_30000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_30000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa


	# ------ 35000
	cam3_infile_35000 = cam3_filebase + 'CliMT_CAM3_1D_35000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
	ds = netCDF4.Dataset(cam3_infile_35000)
	
	cam3_flux_35000_srf = ds.variables['FSWNETSRF'][:]*1.630705885*0.501472063 + ds.variables['FLWNETSRF'][:] - cam3_flux_330_srf
	cam3_flux_35000_200hpa = ds.variables['FSWNET'][19]*1.630705885*0.501472063 + ds.variables['FLWNET'][19] - cam3_flux_330_200hpa
	cam3_flux_35000_toa = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330_toa




	#--------------------------------------------------------
	# SMART
	#--------------------------------------------------------

	# ------ Read in fluxes 330
	smart_infile = smart_filebase + 'fys_co2_330.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 330
	smart_flux_330 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063


	# ------ Read in fluxes 100
	smart_infile = smart_filebase + 'fys_co2_100.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 100
	smart_flux_100 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------ Read in fluxes 500
	smart_infile = smart_filebase + 'fys_co2_500.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 500
	smart_flux_500 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330

	# ------ Read in fluxes 1000
	smart_infile = smart_filebase + 'fys_co2_1000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 1000
	smart_flux_1000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330

	# ------ Read in fluxes 1500
	smart_infile = smart_filebase + 'fys_co2_1500.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 1500
	smart_flux_1500 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------ Read in fluxes 2000
	smart_infile = smart_filebase + 'fys_co2_2000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 2000
	smart_flux_2000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------ Read in fluxes 5000
	smart_infile = smart_filebase + 'fys_co2_5000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 5000
	smart_flux_5000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------ Read in fluxes 10000
	smart_infile = smart_filebase + 'fys_co2_10000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 10000
	smart_flux_10000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------ Read in fluxes 20000
	smart_infile = smart_filebase + 'fys_co2_20000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 20000
	smart_flux_20000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330

	# ------ Read in fluxes 30000
	smart_infile = smart_filebase + 'fys_co2_30000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)	
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])	

	# calc total fluxes 30000
	smart_flux_30000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------ Read in fluxes 35000
	smart_infile = smart_filebase + 'fys_co2_35000.txt'
	dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

	hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)
	sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

	# calc total fluxes 35000
	smart_flux_35000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330


	# ------------------ Organize data for plotting -----------

	# CAM5
	fluxes_cam5_srf = np.array([np.reshape(c5_total_flux_100, 42)[0], np.reshape(c5_total_flux_330, 42)[0] - np.reshape(c5_total_flux_330, 42)[0], np.reshape(c5_total_flux_500, 42)[0], np.reshape(c5_total_flux_1000, 42)[0],np.reshape(c5_total_flux_1500, 42)[0], np.reshape(c5_total_flux_2000, 42)[0], np.reshape(c5_total_flux_5000, 42)[0], np.reshape(c5_total_flux_10000, 42)[0],np.reshape(c5_total_flux_20000, 42)[0], np.reshape(c5_total_flux_30000, 42)[0], np.reshape(c5_total_flux_35000, 42)[0]])
	fluxes_cam5_200hpa = np.array([np.reshape(c5_total_flux_100, 42)[23], np.reshape(c5_total_flux_330, 42)[23] - np.reshape(c5_total_flux_330, 42)[23], np.reshape(c5_total_flux_500, 42)[23], np.reshape(c5_total_flux_1000, 42)[23],np.reshape(c5_total_flux_1500, 42)[23], np.reshape(c5_total_flux_2000, 42)[23], np.reshape(c5_total_flux_5000, 42)[23], np.reshape(c5_total_flux_10000, 42)[23],np.reshape(c5_total_flux_20000, 42)[23], np.reshape(c5_total_flux_30000, 42)[23], np.reshape(c5_total_flux_35000, 42)[23]])
	fluxes_cam5_toa = np.array([np.reshape(c5_total_flux_100, 42)[41], np.reshape(c5_total_flux_330, 42)[41] - np.reshape(c5_total_flux_330, 42)[41], np.reshape(c5_total_flux_500, 42)[41], np.reshape(c5_total_flux_1000, 42)[41],np.reshape(c5_total_flux_1500, 42)[41], np.reshape(c5_total_flux_2000, 42)[41], np.reshape(c5_total_flux_5000, 42)[41], np.reshape(c5_total_flux_10000, 42)[41],np.reshape(c5_total_flux_20000, 42)[41], np.reshape(c5_total_flux_30000, 42)[41], np.reshape(c5_total_flux_35000, 42)[41]])


	# CAM3
	fluxes_cam3_srf = np.array([cam3_flux_100_srf, cam3_flux_330_srf - cam3_flux_330_srf, cam3_flux_500_srf, cam3_flux_1000_srf, cam3_flux_2000_srf, cam3_flux_5000_srf, cam3_flux_10000_srf, cam3_flux_20000_srf, cam3_flux_30000_srf, cam3_flux_35000_srf])
	fluxes_cam3_200hpa = np.array([cam3_flux_100_200hpa, cam3_flux_330_200hpa - cam3_flux_330_200hpa, cam3_flux_500_200hpa, cam3_flux_1000_200hpa, cam3_flux_2000_200hpa, cam3_flux_5000_200hpa, cam3_flux_10000_200hpa, cam3_flux_20000_200hpa, cam3_flux_30000_200hpa, cam3_flux_35000_200hpa])
	fluxes_cam3_toa = np.array([cam3_flux_100_toa, cam3_flux_330_toa - cam3_flux_330_toa, cam3_flux_500_toa, cam3_flux_1000_toa, cam3_flux_2000_toa, cam3_flux_5000_toa, cam3_flux_10000_toa, cam3_flux_20000_toa, cam3_flux_30000_toa, cam3_flux_35000_toa])

	
	# SMART
	fluxes_smart_srf = np.array([smart_flux_100[39], smart_flux_330[39] - smart_flux_330[39], smart_flux_500[39], smart_flux_1000[39], smart_flux_1500[39], smart_flux_2000[39], smart_flux_5000[39], smart_flux_10000[39], smart_flux_20000[39], smart_flux_30000[39], smart_flux_35000[39]])
	fluxes_smart_200hpa = np.array([smart_flux_100[10], smart_flux_330[10] - smart_flux_330[10], smart_flux_500[10], smart_flux_1000[10], smart_flux_1500[10], smart_flux_2000[10], smart_flux_5000[10], smart_flux_10000[10], smart_flux_20000[10], smart_flux_30000[10], smart_flux_35000[10]])
	fluxes_smart_toa = np.array([smart_flux_100[0], smart_flux_330[0] - smart_flux_330[0], smart_flux_500[0], smart_flux_1000[0], smart_flux_1500[0], smart_flux_2000[0], smart_flux_5000[0], smart_flux_10000[0], smart_flux_20000[0], smart_flux_30000[0], smart_flux_35000[0]])


	fluxes_all = [fluxes_cam5_srf, fluxes_cam5_200hpa, fluxes_cam5_toa, 
					fluxes_cam3_srf, fluxes_cam3_200hpa, fluxes_cam3_toa, 
					fluxes_smart_srf, fluxes_smart_200hpa, fluxes_smart_toa]
	
	return(fluxes_all)


# -------------------------------------------------------------------------------
# This function preps the model output for Extended Data figure 6.
# -------------------------------------------------------------------------------
def prep_anomaly_histograms(filebase, filebase_c5, outloc):
	cases_cam4 = {'08','09','10'} # we need the S/S0=0.9 case for these plots
	cases_cam5 = {'09','10'}

	# --- subset the model output for any maps that we don't have yet

	# prep the eis, lts, waterpath maps if haven't yet got it
	map_prep_eis(filebase, outloc, 'cam4')
	map_prep_eis(filebase_c5, outloc, 'cam5')
	map_prep_lts(filebase, outloc, 'cam4')
	map_prep_lts(filebase_c5, outloc, 'cam5')
	map_total_waterpath(filebase, outloc, 'cam4')

	c4_eis_1_file = ''
	c4_eis_8_file = ''
	c4_eis_9_file = ''
	c5_eis_1_file = ''
	c5_eis_9_file = ''

	c4_wp_1_file = ''
	c4_wp_8_file = ''
	c4_wp_9_file = ''

	c4_swcf_1_file = ''
	c4_swcf_8_file = ''
	c4_swcf_9_file = ''
	c5_swcf_1_file = ''
	c5_swcf_9_file = ''

	c4_lhflx_1_file = ''
	c4_lhflx_8_file = ''
	c4_lhflx_9_file = ''
	c5_lhflx_1_file = ''
	c5_lhflx_9_file = ''

	c4_cldlow_1_file = ''
	c4_cldlow_8_file = ''
	c4_cldlow_9_file = ''
	c5_cldlow_1_file = ''
	c5_cldlow_9_file = ''

	c4_ah_1_file = ''
	c4_ah700_1_file = ''
	c4_ah1000_1_file = ''

	c4_ah_8_file = ''
	c4_ah700_8_file = ''
	c4_ah1000_8_file = ''
	
	c4_ah_9_file = ''
	c4_ah700_9_file = ''
	c4_ah1000_9_file = ''
	c5_ah_1_file = ''
	c5_ah700_1_file = ''
	c5_ah1000_1_file = ''
	
	c5_ah_9_file = ''
	c5_ah700_9_file = ''
	c5_ah1000_9_file = ''

	c4_ts_1_file = ''
	c4_ts_8_file = ''
	c4_ts_9_file = ''
	c5_ts_1_file = ''
	c5_ts_9_file = ''
	
	# ------ CAM4
	cam_version = 'cam4'
	for case in cases_cam4:
		selyear = '-selyear,21/40'	
		sellevel = '-sellevel,696.79629'

		eis = outloc + 'c4_eis_' + case + '.nc'
		eis = outloc + 'c4_waterpath_' + case + '.nc'

		# for absolute humidity
		ah700 = outloc + 'c4_700ah_' + case + '.nc'
		ah1000 = outloc  + 'c4_1000ah_' + case + '.nc'
	
		# for others
		swcf = outloc + 'c4_swcf_' + case + '.nc'
		cldlow = outloc + 'c4_cldlow_' + case + '.nc'
		lhflx = outloc + 'c4_lhflx_' + case + '.nc'

		# for ts
		ts = outloc + 'c4_ts_' + case + '.nc'

		if case == '08':
			c4_eis_8_file = eis
			c4_wp_8_file = wp
			c4_swcf_8_file = swcf
			c4_cldlow_8_file
			c4_ah700_8_file = ah700
			c4_ah1000_8_file = ah700
			c4_ts_8_file = ts


		# absolute humidity
		if not os.path.isfile(ah700):
			infile = filebase + cam_version + '_' + case + '.nc'

			# calc Q at 700 hPa
			syscall = r"//usr//bin//cdo timmean " + selyear + ' -select,name=Q ' + sellevel +  ' ' + infile + ' ' + ah700
			print(syscall)
			os.system(syscall)

		# absolute humidity
		if not os.path.isfile(ah1000):
			infile = filebase + cam_version + '_' + case + '.nc'

			# calc Q at surface
			syscall = r"//usr//bin//cdo timmean " + selyear + ' -select,name=QREFHT ' + infile + ' ' + ah1000
			print(syscall)
			os.system(syscall)

		# Other variables
		if not os.path.isfile(swcf):
			infile = filebase + cam_version + '_' + case + '.nc'
			# get swcf
			syscall = r"//usr//bin//cdo  timmean " + selyear + ' -select,name=SWCF '+infile + ' ' + swcf
			print(syscall)
			os.system(syscall)



		# CLDLOW
		if not os.path.isfile(cldlow):
			infile = filebase + cam_version + '_' + case + '.nc'

			# get cloudlow
			syscall = r"//usr//bin//cdo  timmean "  + selyear + ' -select,name=CLDLOW '+infile + ' ' + cldlow
			print(syscall)
			os.system(syscall)


		# LHFLX
		if not os.path.isfile(lhflx):
			infile = filebase + cam_version + '_' + case + '.nc'
			
			# get lhflx
			syscall = r"//usr//bin//cdo  timmean " + selyear + ' -select,name=LHFLX '+infile + ' ' + lhflx
			print(syscall)
			os.system(syscall)


		# TS
		if not os.path.isfile(ts):
			infile = filebase + cam_version + '_' + case + '.nc'

			# get ts
			syscall = r"//usr//bin//cdo  timmean " + selyear + ' -select,name=TS '+infile + ' ' + ts
			print(syscall)
			os.system(syscall)


	# ----- CAM5
	cam_version='cam5'
	for case in cases_cam5:
		selyear = '-selyear,31/40'	
		sellevel = '-sellevel,691.389430314302'


		# for absolute humidity
		ah700 = outloc + 'c5_700ah_' + case + '.nc'
		ah1000 = outloc  + 'c5_1000ah_' + case + '.nc'
	
		# for others
		swcf = outloc + 'c5_swcf_' + case + '.nc'
		cldlow = outloc + 'c5_cldlow_' + '.nc'
		lhflx = outloc + 'c5_lhflx_' + '.nc'

		# for ts
		ts = outloc + 'c5_ts_' + case + '.nc'


		# absolute humidity
		if not os.path.isfile(ah700):
			infile = filebase_c5 + cam_version + '_' + case + '.nc'

			# calc Q at 700 hPa
			syscall = r"//usr//bin//cdo timmean " + selyear + ' -select,name=Q ' + sellevel +  ' ' + infile + ' ' + ah700
			print(syscall)
			os.system(syscall)
		
		# absolute humidity
		if not os.path.isfile(ah1000):
			infile = filebase_c5 + cam_version + '_' + case + '.nc'

			# calc Q at surface
			syscall = r"//usr//bin//cdo timmean " + selyear + ' -select,name=QREFHT ' + infile + ' ' + ah1000
			print(syscall)
			os.system(syscall)


		# Other variables
		if not os.path.isfile(swcf):
			infile = filebase_c5 + cam_version + '_' + case + '.nc'
			# get swcf
			syscall = r"//usr//bin//cdo  timmean " + selyear + ' -select,name=SWCF '+infile + ' ' + swcf
			print(syscall)
			os.system(syscall)



		# CLDLOW
		if not os.path.isfile(cldlow):
			infile = filebase_c5 + cam_version + '_' + case + '.nc'

			# get cloudlow
			syscall = r"//usr//bin//cdo  timmean "  + selyear + ' -select,name=CLDLOW '+infile + ' ' + cldlow
			print(syscall)
			os.system(syscall)


		# LHFLX
		if not os.path.isfile(lhflx):
			infile = filebase_c5 + cam_version + '_' + case + '.nc'
			
			# get lhflx
			syscall = r"//usr//bin//cdo  timmean " + selyear + ' -select,name=LHFLX '+infile + ' ' + lhflx
			print(syscall)
			os.system(syscall)


		# TS
		if not os.path.isfile(ts):
			infile = filebase_c5 + cam_version + '_' + case + '.nc'

			# get ts
			syscall = r"//usr//bin//cdo  timmean " + selyear + ' -select,name=TS '+infile + ' ' + ts
			print(syscall)
			os.system(syscall)


	# --- calculate the anomalies and put the data in new .nc files

	