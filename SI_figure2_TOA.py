#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import climt
import netCDF4
import cdo
import json

download_path = '' # enter the path to the directory where you downloaded the archived data, eg '/home/user/Downloads'

cam5_filebase = download_path + '/FYSP_clouds_archive/radiative-transfer/CAM5/'
cam3_filebase = download_path + '/FYSP_clouds_archive/radiative-transfer/CAM3/'
smart_filebase = download_path + '/FYSP_clouds_archive/radiative-transfer/SMART/'

# -------------------------     CAM5     -----------------------------
# Read in fluxes
cam5_infile = cam5_filebase + 'fluxes_100.txt'
ds = json.load(open(cam5_infile))

c5lwup_100 = np.asarray(ds['lw_up'])
c5lwdn_100 = np.asarray(ds['lw_dn'])
c5swup_100 = np.asarray(ds['sw_up'])
c5swdn_100 = np.asarray(ds['sw_dn'])
#net_100 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_330.txt'
ds = json.load(open(cam5_infile))

c5lwup_330 = np.asarray(ds['lw_up'])
c5lwdn_330 = np.asarray(ds['lw_dn'])
c5swup_330 = np.asarray(ds['sw_up'])
c5swdn_330 = np.asarray(ds['sw_dn'])
#net_330 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_500.txt'
ds = json.load(open(cam5_infile))

c5lwup_500 = np.asarray(ds['lw_up'])
c5lwdn_500 = np.asarray(ds['lw_dn'])
c5swup_500 = np.asarray(ds['sw_up'])
c5swdn_500 = np.asarray(ds['sw_dn'])
#net_500 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_1000.txt'
ds = json.load(open(cam5_infile))

c5lwup_1000 = np.asarray(ds['lw_up'])
c5lwdn_1000 = np.asarray(ds['lw_dn'])
c5swup_1000 = np.asarray(ds['sw_up'])
c5swdn_1000 = np.asarray(ds['sw_dn'])
#net_1000 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_1500.txt'
ds = json.load(open(cam5_infile))

c5lwup_1500 = np.asarray(ds['lw_up'])
c5lwdn_1500 = np.asarray(ds['lw_dn'])
c5swup_1500 = np.asarray(ds['sw_up'])
c5swdn_1500 = np.asarray(ds['sw_dn'])
#net_1500 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_2000.txt'
ds = json.load(open(cam5_infile))

c5lwup_2000 = np.asarray(ds['lw_up'])
c5lwdn_2000 = np.asarray(ds['lw_dn'])
c5swup_2000 = np.asarray(ds['sw_up'])
c5swdn_2000 = np.asarray(ds['sw_dn'])
#net_2000 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_5000.txt'
ds = json.load(open(cam5_infile))

c5lwup_5000 = np.asarray(ds['lw_up'])
c5lwdn_5000 = np.asarray(ds['lw_dn'])
c5swup_5000 = np.asarray(ds['sw_up'])
c5swdn_5000 = np.asarray(ds['sw_dn'])
#net_5000 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_10000.txt'
ds = json.load(open(cam5_infile))

c5lwup_10000 = np.asarray(ds['lw_up'])
c5lwdn_10000 = np.asarray(ds['lw_dn'])
c5swup_10000 = np.asarray(ds['sw_up'])
c5swdn_10000 = np.asarray(ds['sw_dn'])
#net_10000 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_20000.txt'
ds = json.load(open(cam5_infile))

c5lwup_20000 = np.asarray(ds['lw_up'])
c5lwdn_20000 = np.asarray(ds['lw_dn'])
c5swup_20000 = np.asarray(ds['sw_up'])
c5swdn_20000 = np.asarray(ds['sw_dn'])
#net_20000 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_30000.txt'
ds = json.load(open(cam5_infile))

c5lwup_30000 = np.asarray(ds['lw_up'])
c5lwdn_30000 = np.asarray(ds['lw_dn'])
c5swup_30000 = np.asarray(ds['sw_up'])
c5swdn_30000 = np.asarray(ds['sw_dn'])
#net_30000 = np.asarray(ds['net'])

cam5_infile = cam5_filebase + 'fluxes_35000.txt'
ds = json.load(open(cam5_infile))

c5lwup_35000 = np.asarray(ds['lw_up'])
c5lwdn_35000 = np.asarray(ds['lw_dn'])
c5swup_35000 = np.asarray(ds['sw_up'])
c5swdn_35000 = np.asarray(ds['sw_dn'])
#net_35000 = np.asarray(ds['net'])


# ------------ Calc total fluxes CAM5 --------------
total_flux_330 = c5swup_330*(-1) + c5swdn_330 - c5lwup_330 + c5lwdn_330
total_flux_100 = c5swup_100*(-1) + c5swdn_100 - c5lwup_100 + c5lwdn_100 - total_flux_330
total_flux_500 = c5swup_500*(-1) + c5swdn_500 - c5lwup_500 + c5lwdn_500 - total_flux_330
total_flux_1000 = c5swup_1000*(-1) + c5swdn_1000 - c5lwup_1000 + c5lwdn_1000 - total_flux_330
total_flux_1500 = c5swup_1500*(-1) + c5swdn_1500 - c5lwup_1500 + c5lwdn_1500- total_flux_330
total_flux_2000 = c5swup_2000*(-1) + c5swdn_2000 - c5lwup_2000 + c5lwdn_2000 - total_flux_330
total_flux_5000 = c5swup_5000*(-1) + c5swdn_5000 - c5lwup_5000 + c5lwdn_5000 - total_flux_330
total_flux_10000 = c5swup_10000*(-1) + c5swdn_10000 - c5lwup_10000 + c5lwdn_10000 - total_flux_330
total_flux_20000 = c5swup_20000*(-1) + c5swdn_20000 - c5lwup_20000 + c5lwdn_20000 - total_flux_330
total_flux_30000 = c5swup_30000*(-1) + c5swdn_30000 - c5lwup_30000 + c5lwdn_30000 - total_flux_330
total_flux_35000 = c5swup_35000*(-1) + c5swdn_35000 - c5lwup_35000 + c5lwdn_35000 - total_flux_330




#------------------------ CAM3 --------------------------------
# 330
cam3_infile_330 = cam3_filebase + '/CliMT_CAM3_1D_330.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_330)
cam3_flux_330 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:]

cam3_pressure = ds.variables['Pres'][:]

# 100
cam3_infile_100 = cam3_filebase + 'CliMT_CAM3_1D_100.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_100)
cam3_flux_100 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 500
cam3_infile_500 = cam3_filebase + 'CliMT_CAM3_1D_500.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_500)
cam3_flux_500 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 1000
cam3_infile_1000 = cam3_filebase + 'CliMT_CAM3_1D_1000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_1000)
cam3_flux_1000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330


# 2000
cam3_infile_2000 = cam3_filebase + 'CliMT_CAM3_1D_2000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_2000)
cam3_flux_2000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 5000
cam3_infile_5000 = cam3_filebase + 'CliMT_CAM3_1D_5000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_5000)
cam3_flux_5000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 10000
cam3_infile_10000 = cam3_filebase + 'CliMT_CAM3_1D_10000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_10000)
cam3_flux_10000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 20000
cam3_infile_20000 = cam3_filebase + 'CliMT_CAM3_1D_20000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_20000)
cam3_flux_20000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 30000
cam3_infile_30000 = cam3_filebase + 'CliMT_CAM3_1D_30000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_30000)
cam3_flux_30000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330

# 35000
cam3_infile_35000 = cam3_filebase + 'CliMT_CAM3_1D_35000.0ppmvCO2_1842ppbvCH4_328ppbvN2O_42levels.nc'
ds = netCDF4.Dataset(cam3_infile_35000)
cam3_flux_35000 = ds.variables['FSWNETTOA'][:]*1.630705885*0.501472063 + ds.variables['FLWNETTOA'][:] - cam3_flux_330





# ----------------------- SMART -----------------------------


# Read in fluxes 330
smart_infile = smart_filebase + 'fys_co2_330.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 330
smart_flux_330 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063



# Read in fluxes 100
smart_infile = smart_filebase + 'fys_co2_100.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 100
smart_flux_100 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330



# Read in fluxes 500
smart_infile = smart_filebase + 'fys_co2_500.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 500
smart_flux_500 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330




# Read in fluxes 1000
smart_infile = smart_filebase + 'fys_co2_1000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 1000
smart_flux_1000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330




# Read in fluxes 1500
smart_infile = smart_filebase + 'fys_co2_1500.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 1500
smart_flux_1500 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330





# Read in fluxes 2000
smart_infile = smart_filebase + 'fys_co2_2000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 2000
smart_flux_2000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330



# Read in fluxes 5000
smart_infile = smart_filebase + 'fys_co2_5000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 5000
smart_flux_5000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330



# Read in fluxes 10000
smart_infile = smart_filebase + 'fys_co2_10000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 10000
smart_flux_10000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330




# Read in fluxes 20000
smart_infile = smart_filebase + 'fys_co2_20000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 20000
smart_flux_20000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330




# Read in fluxes 30000
smart_infile = smart_filebase + 'fys_co2_30000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 30000
smart_flux_30000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330



# Read in fluxes 35000
smart_infile = smart_filebase + 'fys_co2_35000.txt'
dtype2 = np.dtype([('Pressure_hrt','f8'),('Temperature','f8'),('Altitude','f8'),('solarq','f8'),('thermalq','f8'),('dirsolarflux','f8'),('dnsolarflx','f8'),('upsolarflx','f8'),('dnthermalflx','f8'),('upthermalflx','f8'),('pflux','f8'),('tflux','f8')])

hrt_profile = np.loadtxt(smart_infile, dtype=dtype2, skiprows=1)

sw_down_hrt = np.array(hrt_profile['dnsolarflx']) + np.array(hrt_profile['dirsolarflux'])

# calc total fluxes 35000
smart_flux_35000 = np.array(hrt_profile['dnthermalflx']) - np.array(hrt_profile['upthermalflx']) + sw_down_hrt*0.501472063 - np.array(hrt_profile['upsolarflx'])*0.501472063 - smart_flux_330



# ------------------ Organize data for plotting -----------

# CAM5
fluxes_cam5 = np.array([np.reshape(total_flux_100, 42)[41], np.reshape(total_flux_330, 42)[41] - np.reshape(total_flux_330, 42)[41], np.reshape(total_flux_500, 42)[41], np.reshape(total_flux_1000, 42)[41],np.reshape(total_flux_1500, 42)[41], np.reshape(total_flux_2000, 42)[41], np.reshape(total_flux_5000, 42)[41], np.reshape(total_flux_10000, 42)[41],np.reshape(total_flux_20000, 42)[41], np.reshape(total_flux_30000, 42)[41], np.reshape(total_flux_35000, 42)[41]])

#fluxes_cam5 = np.array([net_100[41], net_330[41] - net_330[41], net_500[41], net_1000[41], net_1500[41], net_2000[41], net_5000[41], net_10000[41], net_20000[41], net_30000[41], net_35000[41]])



# CAM3
fluxes_cam3 = np.array([cam3_flux_100, cam3_flux_330 - cam3_flux_330, cam3_flux_500, cam3_flux_1000, cam3_flux_2000, cam3_flux_5000, cam3_flux_10000, cam3_flux_20000, cam3_flux_30000, cam3_flux_35000])

# SMART
fluxes_smart = np.array([smart_flux_100[0], smart_flux_330[0] - smart_flux_330[0], smart_flux_500[0], smart_flux_1000[0], smart_flux_1500[0], smart_flux_2000[0], smart_flux_5000[0], smart_flux_10000[0], smart_flux_20000[0], smart_flux_30000[0], smart_flux_35000[0]])


co2_cam5 = np.array([100, 330, 500, 1000, 1500, 2000, 5000, 10000, 20000, 30000, 35000])
co2_cam3 = np.array([100, 330, 500, 1000, 2000, 5000, 10000, 20000, 30000, 35000])
co2_smart = np.array([100, 330, 500, 1000, 1500, 2000, 5000, 10000, 20000, 30000, 35000])

# Make plots
f = plt.figure()

plt.plot(co2_cam5, fluxes_cam5, color='blue', label='CAM5')
plt.plot(co2_cam3, fluxes_cam3, color='green', label='CAM3')
plt.plot(co2_smart, fluxes_smart, color='black', label='SMART')
plt.title('Radiative Forcing - Top of Atmosphere')
plt.xlabel('CO2')
plt.ylabel('Forcing (W/m2)')

# Display the legend below the plot
plt.legend()

plt.xscale("log")

plt.show()
f.savefig('rad_forcing_TOA_log.pdf', bbox_inches='tight')
