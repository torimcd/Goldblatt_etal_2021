#!/usr/bin/env python3
"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""

import os
import sys
import numpy
import netCDF4
import operator
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

run = 'Final'
years = '40'

sw_dn = 'FSDS'
sw_net = 'FSNS'
lw_dn = 'FLDS'
lw_net = 'FLNS'

lhflx = 'LHFLX'
shflx = 'SHFLX'

filebase = '/home/vmcd/projects/modelOutput/gcmexpt/CAM4/'
filebase_c5 = '/home/vmcd/projects/modelOutput/gcmexpt/CAM5/'
outfilebase = 'srfenergy_'

casenames = ['1.1','1.075','1.05','1.025','1.0','0.975','0.95','0.925','0.9','0.875','0.85','0.825','0.8','0.775','0.75','0.725','0.7']
casenames_c5 = ['1.05','1.025','1.0', '0.975', '0.95','0.925','0.9']

sc_all = [1.1, 1.075, 1.05, 1.025, 1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.725, 0.7]
sc_c5 = ['1.05','1.025','1.0', '0.975', '0.95','0.925','0.9']

for CASENAME in casenames:
	outfile_case = filebase+CASENAME+'/'+outfilebase+CASENAME+'.nc'
	# check directly if the file exists
	if not os.path.isfile(outfile_case):
		if os.path.isdir(filebase+CASENAME):
			infile = filebase + CASENAME+'/merged.nc'
			# calc cldlow global average per month
			syscall = 'cdo timmean -fldmean -selyear,21/40 -select,name=FSDS,FSDSC,FSNS,FSNSC,FLDS,FLDSC,FLNS,FLNSC,LHFLX,SHFLX '+infile+ ' ' +outfile_case
			os.system(syscall)

for CASENAME in casenames_c5:
	outfile_case = filebase_c5+CASENAME+'/'+outfilebase+CASENAME+'.nc'
	# check directly if the file exists
	if not os.path.isfile(outfile_case):
		if os.path.isdir(filebase_c5+CASENAME):
			infile = filebase_c5 + CASENAME+'/merged.nc'
			# calc cldlow global average per month
			syscall = 'cdo timmean -fldmean -selyear,31/60 -select,name=FSDS,FSDSC,FSNS,FSNSC,FLDS,FLDSC,FLNS,FLNSC,LHFLX,SHFLX '+infile+ ' ' +outfile_case
			os.system(syscall)


#create plot
fig = plt.figure()
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
	dsloc = filebase+CASENAME+'/'+outfilebase+CASENAME+'.nc'
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
	dsloc = filebase_c5+CASENAME+'/'+outfilebase+CASENAME+'.nc'
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
plt.legend(handles2, labels2, bbox_to_anchor=(0., -0.3, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.)

plt.title('Global Average Surface Energy Budget', fontsize=16)
plt.xlabel(r'$\mathsf{S/S_0}$', fontsize=16)
plt.ylabel(r'$\mathsf{Flux}$' + ' ' +r'$\mathsf{(W/m^2)}$', fontsize=16)
plt.minorticks_on()
plt.axis([0.675,1.125,0,400])
plt.grid()
plt.show()

fig.savefig("srf_energy_budget.pdf", bbox_inches='tight')

