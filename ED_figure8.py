#!/usr/bin/env python3

"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

This script plots the radiative forcing for all models + SMART.

"""
import matplotlib
matplotlib.use("Agg")

import os
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import processing_functions as pf


# ------------------------------------------------------------------------
# change this section to match where you downloaded the model output files 
# ------------------------------------------------------------------------

download_path = '' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'
filebase = download_path + 'FYSP_clouds_archive/radiative_transfer/'


# ----------------------------------------------------
# calculate the fluxes for all models for all levels
# ----------------------------------------------------
fluxes = pf.calculate_rad_fluxes(filebase)

cam5_fluxes_srf = fluxes[0]
cam5_fluxes_200hpa = fluxes[1]
cam5_fluxes_toa = fluxes[2]

cam3_fluxes_srf = fluxes[3]
cam3_fluxes_200hpa = fluxes[4]
cam3_fluxes_toa = fluxes[5]

smart_fluxes_srf = fluxes[6]
smart_fluxes_200hpa = fluxes[7]
smart_fluxes_toa = fluxes[8]


co2_cam5 = np.array([100, 330, 500, 1000, 1500, 2000, 5000, 10000, 20000, 30000, 35000])
co2_cam3 = np.array([100, 330, 500, 1000, 2000, 5000, 10000, 20000, 30000, 35000])
co2_smart = np.array([100, 330, 500, 1000, 1500, 2000, 5000, 10000, 20000, 30000, 35000])



# ------------------
# make the figure
# ------------------
fig = plt.figure(figsize=(3.46457, 5))
#fig = plt.figure(figsize=(8.5,11))
ax1 = fig.add_subplot(311)

ax1.plot(co2_cam5, cam5_fluxes_toa, color='blue', label='CAM5')
ax1.plot(co2_cam3, cam3_fluxes_toa, color='green', label='CAM3')
ax1.plot(co2_smart, smart_fluxes_toa, color='black', label='SMART')
ax1.set_title('Radiative Forcing - Top of Atmosphere', fontsize=7)
ax1.set_xlabel(r'$\mathsf{CO_2}$' + ' (ppmv)', fontsize=7)
ax1.set_ylabel('Forcing ' + r'$\mathsf{W/m^2}$', fontsize=7)
ax1.set_xscale("log")
ax1.tick_params(labelsize=6)
plt.legend(loc='lower right', prop={"size":5})

ax2 = fig.add_subplot(312)

ax2.plot(co2_cam5, cam5_fluxes_200hpa, color='blue', label='CAM5')
ax2.plot(co2_cam3, cam3_fluxes_200hpa, color='green', label='CAM3')
ax2.plot(co2_smart, smart_fluxes_200hpa, color='black', label='SMART')
ax2.set_title('Radiative Forcing - 200hpa', fontsize=7)
ax2.set_xlabel(r'$\mathsf{CO_2}$' + ' (ppmv)', fontsize=7)
ax2.set_ylabel('Forcing ' + r'$\mathsf{W/m^2}$', fontsize=7)
ax2.set_xscale("log")
ax2.tick_params(labelsize=5)
plt.legend(loc='lower right', prop={"size":5})

ax3 = fig.add_subplot(313)

ax3.plot(co2_cam5, cam5_fluxes_srf, color='blue', label='CAM5')
ax3.plot(co2_cam3, cam3_fluxes_srf, color='green', label='CAM3')
ax3.plot(co2_smart, smart_fluxes_srf, color='black', label='SMART')
ax3.set_title('Radiative Forcing - Surface', fontsize=7)
ax3.set_xlabel(r'$\mathsf{CO_2}$' + ' (ppmv)', fontsize=7)
ax3.set_ylabel('Forcing ' + r'$\mathsf{W/m^2}$', fontsize=7)
ax3.set_xscale("log")
ax3.tick_params(labelsize=5)
plt.legend(loc='lower right', prop={"size":5})

plt.tight_layout()

plt.show()

fig.savefig("figures_ED/ED_figure8.pdf", format='pdf', bbox_inches='tight')

