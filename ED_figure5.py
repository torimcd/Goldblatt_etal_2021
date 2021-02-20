#!/usr/bin/env python3

"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

This script plots the surface energy balance for both CAM4 and CAM5 model output.

"""
import matplotlib
matplotlib.use("Agg")

import os
import sys
import numpy
import netCDF4
import operator
import matplotlib
import matplotlib.pyplot as plt
import processing_functions as pf


# ------------------------------------------------------------------------
# change this section to match where you downloaded the model output files 
# ------------------------------------------------------------------------

download_path = '' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

filebase = download_path + 'FYSP_clouds_archive/CAM4/'
filebase_c5 = download_path + 'FYSP_clouds_archive/CAM5/'
outfileloc = download_path + 'temp_data/' # this is the location to save the processed netcdf files to

# ------------------------------------------


# field variables
fields = 'co2vmr,LHFLX,SHFLX,CLDLOW,CLDHGH,LWCF,SWCF,FSDS,FSDSC,FSNS,FSNSC,FLDS,FLDSC,FLNS,FLNSC'

pf.global_annual_average(filebase, outfileloc, fields, 'cam4')
pf.global_annual_average(filebase_c5, outfileloc, fields, 'cam5')

sw_dn = 'FSDS'
sw_net = 'FSNS'
lw_dn = 'FLDS'
lw_net = 'FLNS'

lhflx = 'LHFLX'
shflx = 'SHFLX'


outfilebase_c4 = 'c4_global_average_'
outfilebase_c5 = 'c5_global_average_'

casenames = ['07','0725','075','0775','08','0825','085','0875','09','0925','095','0975','10','1025','105','1075','11']
casenames_c5 = ['09','0925','095','0975','10','1025','105']

sc_all = [1.1, 1.075, 1.05, 1.025, 1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.725, 0.7]
sc_c5 = ['1.05','1.025','1.0', '0.975', '0.95','0.925','0.9']



#create plot
fig = plt.figure(figsize=(7.08661, 7.28346))
#fig = plt.figure(figsize=(10, 11))
ax = plt.subplot(211)

sw_dn_plot = []
sw_up_plot = []
sw_net_plot = []
lw_dn_plot = []
lw_up_plot = []
lw_net_plot = []

lhflx_plot = []
shflx_plot = []

sw_dn_plot_c5 = []
sw_up_plot_c5 = []
sw_net_plot_c5 = []
lw_dn_plot_c5 = []
lw_up_plot_c5 = []
lw_net_plot_c5 = []


lhflx_plot_c5 = []
shflx_plot_c5 = []

i=0
#plot the data
for CASE in casenames:
	CASENAME = casenames[i]
	dsloc = outfileloc + outfilebase_c4 + CASENAME + '.nc'
	if os.path.isfile(dsloc):
		swd = []
		swdc = []
		swn = []
		swnc = []
		lwd = []
		lwdc = []
		lwn = []
		lwnc = []
		lh = []
		sh = []

		# open the merged file and get out some variables
		dsyear = netCDF4.Dataset(dsloc)
		swd = dsyear.variables[sw_dn][:]
		swn = dsyear.variables[sw_net][:]
		lwd = dsyear.variables[lw_dn][:]
		lwn = dsyear.variables[lw_net][:]
		sh = dsyear.variables[shflx][:]
		lh = dsyear.variables[lhflx][:]
		dsyear.close() #close the file

		swd = swd.flatten()
		swn = swn.flatten()
		lwd = lwd.flatten()
		lwn = lwn.flatten()
		sh = sh.flatten()
		lh = lh.flatten()

		
	sw_dn_plot.append(swd.item(0))
	sw_net_plot.append(swn.item(0))
	sw_up_plot.append((swn.item(0) - swd.item(0))*(-1))
	lw_dn_plot.append(lwd.item(0))
	lw_net_plot.append(lwn.item(0))
	lw_up_plot.append((lwn.item(0) - lwd.item(0))*(-1))
	lhflx_plot.append(lh.item(0))
	shflx_plot.append(sh.item(0))
	i+=1

i=0
for CASENAME in casenames_c5:
	CASENAME = casenames_c5[i]
	dsloc = outfileloc + outfilebase_c5 + CASENAME + '.nc'
	if os.path.isfile(dsloc):
		swd_c5 = []
		swdc_c5 = []
		swn_c5 = []
		swnc_c5 = []
		lwd_c5 = []
		lwdc_c5 = []
		lwn_c5 = []
		lwnc_c5 = []
		lh_c5 = []
		sh_c5 = []

		# open the merged file and get out some variables
		dsyear = netCDF4.Dataset(dsloc)
		swd_c5 = dsyear.variables[sw_dn][:]
		swn_c5 = dsyear.variables[sw_net][:]
		lwd_c5 = dsyear.variables[lw_dn][:]
		lwn_c5 = dsyear.variables[lw_net][:]
		sh_c5 = dsyear.variables[shflx][:]
		lh_c5 = dsyear.variables[lhflx][:]
		dsyear.close() #close the file

		swd_c5 = swd_c5.flatten()
		swn_c5 = swn_c5.flatten()
		lwd_c5 = lwd_c5.flatten()
		lwn_c5 = lwn_c5.flatten()
		sh_c5 = sh_c5.flatten()
		lh_c5 = lh_c5.flatten()

		
	sw_dn_plot_c5.append(swd_c5.item(0))
	sw_net_plot_c5.append(swn_c5.item(0))
	sw_up_plot_c5.append((swn_c5.item(0) - swd_c5.item(0))*(-1))
	lw_dn_plot_c5.append(lwd_c5.item(0))
	lw_net_plot_c5.append(lwn_c5.item(0))
	lw_up_plot_c5.append((lwn_c5.item(0) - lwd_c5.item(0))*(-1))

	lhflx_plot_c5.append(lh_c5.item(0))
	shflx_plot_c5.append(sh_c5.item(0))
	i+=1



plt.plot(sc_all, sw_dn_plot, color='gold', marker='v', linestyle='-', label='Shortwave Down CAM4')
plt.plot(sc_all, sw_up_plot, color='olive', marker='^', linestyle='-', label='Shortwave Up CAM4')
plt.plot(sc_all, sw_net_plot, color='orange', marker='o', label='Shortwave Net CAM4')
plt.plot(sc_all, lw_dn_plot, color='dodgerblue', marker='v', linestyle='-', label='Longwave Down CAM4')
plt.plot(sc_all, lw_up_plot, color='skyblue', marker='^', linestyle='-', label='Longwave Up CAM4')
plt.plot(sc_all, lw_net_plot, color='blue', marker='o', label='Longwave Net CAM4')
plt.plot(sc_all, lhflx_plot, color='green', marker='s', linestyle='-', label='Latent Heat CAM4')
plt.plot(sc_all, shflx_plot, color='red', marker='s', linestyle='-', alpha=0.7, label='Sensible Heat CAM4')



plt.plot(sc_c5, sw_dn_plot_c5, color='black', marker='x', linestyle='--', label='CAM5')
plt.plot(sc_c5, sw_up_plot_c5, color='black', marker='x', linestyle='--')
plt.plot(sc_c5, sw_net_plot_c5, color='black', marker='x', linestyle='--')
plt.plot(sc_c5, lw_dn_plot_c5, color='black', marker='x', linestyle='--')
plt.plot(sc_c5, lw_up_plot_c5, color='black', marker='x', linestyle='--')
plt.plot(sc_c5, lw_net_plot_c5, color='black', marker='x', linestyle='--')
plt.plot(sc_c5, lhflx_plot_c5, color='black', marker='x', linestyle='--')
plt.plot(sc_c5, shflx_plot_c5, color='black', marker='x', linestyle='--')



#sort legend
handles, labels = ax.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles2, labels2 = zip(*hl)

# Display the legend below the plot
plt.legend(handles2, labels2, bbox_to_anchor=(0., -0.3, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0., prop={"size":7})

plt.title('Global Average Surface Energy Budget', fontsize=7)
plt.xlabel(r'$\mathsf{S/S_0}$', fontsize=7)
plt.ylabel(r'$\mathsf{Flux}$' + ' ' +r'$\mathsf{(W/m^2)}$', fontsize=7)
plt.minorticks_on()
ax.tick_params(labelsize=7) 
plt.axis([0.675,1.125,0,400])
plt.grid()
plt.show()

fig.savefig("figures_ED/ED_figure5.pdf", format='pdf', bbox_inches='tight')

