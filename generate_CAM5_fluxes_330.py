"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""

import numpy as np
import matplotlib.pyplot as plt
import climt
import json

# -------------------------     CAM5     -----------------------------
# Read in pressure/temp profile
infile = 'rad_tr_profile.txt'
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


#------------------------- Shortwave -----------------------------
swradiation = climt.RRTMGShortwave()
surface = climt.SlabSurface()

# Create shortwave state
m_levs=dict(label='levels',values=np.arange(41), units='')
i_levs=dict(label='levels', values=np.arange(42), units='')
swstate=climt.get_default_state([swradiation, surface], mid_levels=m_levs, interface_levels=i_levs)

# Change values to match my atmosphere profile
swstate['surface_temperature'] = swstate['surface_temperature'] - 11.76
swstate['air_pressure'].values = p_mid.reshape(1,1,41)
swstate['air_pressure_on_interface_levels'].values = p_int.reshape(1,1,42)
swstate['air_temperature'].values = T_mid.reshape(1,1,41)
swstate['specific_humidity'].values = Water_mid.reshape(1,1,41)
swstate['mole_fraction_of_ozone_in_air'].values = Ozone_mid.reshape(1,1,41)
swstate['surface_albedo_for_direct_shortwave'] = swstate['surface_albedo_for_direct_shortwave'] +0.06
swstate['single_scattering_albedo_due_to_cloud'] = swstate['single_scattering_albedo_due_to_cloud'] +0.06
swstate['surface_albedo_for_diffuse_shortwave'] = swstate['surface_albedo_for_diffuse_shortwave'] +0.06
swstate['single_scattering_albedo_due_to_aerosol'] = swstate['single_scattering_albedo_due_to_aerosol'] +0.06
swstate['surface_albedo_for_diffuse_near_infrared'] = swstate['surface_albedo_for_diffuse_near_infrared'] +0.06
swstate['surface_albedo_for_direct_near_infrared'] = swstate['surface_albedo_for_direct_near_infrared'] +0.06
swstate['zenith_angle'] = swstate['zenith_angle'] + np.pi/3

swstate['mole_fraction_of_methane_in_air'] = swstate['mole_fraction_of_methane_in_air'] + 0.00000184278
swstate['mole_fraction_of_nitrous_oxide_in_air'] = swstate['mole_fraction_of_nitrous_oxide_in_air'] + 3.28*10**-7

swstate['air_pressure'].attrs['units'] = 'hPa'
swstate['air_pressure_on_interface_levels'].attrs['units'] = 'hPa'

# See the state of all variables
#swstate

# Calc unperturbed fluxes - CO2 at 330 ppm default
swtendencies, swdiagnostics = swradiation(swstate)

unpert_swup = swdiagnostics['upwelling_shortwave_flux_in_air']
unpert_swdown = swdiagnostics['downwelling_shortwave_flux_in_air']


# ----------------- Longwave -------------------------
lwradiation = climt.RRTMGLongwave()
surface = climt.SlabSurface()

# Create shortwave state
m_levs=dict(label='levels',values=np.arange(41), units='')
i_levs=dict(label='levels', values=np.arange(42), units='')
lwstate=climt.get_default_state([lwradiation, surface], mid_levels=m_levs, interface_levels=i_levs)

# Change values to match my atmosphere profile
lwstate['surface_temperature'] = lwstate['surface_temperature'] - 11.76
lwstate['air_pressure'].values = p_mid.reshape(1,1,41)
lwstate['air_pressure_on_interface_levels'].values = p_int.reshape(1,1,42)
lwstate['air_temperature'].values = T_mid.reshape(1,1,41)
lwstate['specific_humidity'].values = Water_mid.reshape(1,1,41)
lwstate['mole_fraction_of_ozone_in_air'].values = Ozone_mid.reshape(1,1,41)

lwstate['mole_fraction_of_methane_in_air'] = lwstate['mole_fraction_of_methane_in_air'] + 0.00000184278
lwstate['mole_fraction_of_nitrous_oxide_in_air'] = lwstate['mole_fraction_of_nitrous_oxide_in_air'] + 3.28*10**-7

lwstate['air_pressure'].attrs['units'] = 'hPa'
lwstate['air_pressure_on_interface_levels'].attrs['units'] = 'hPa'

# See the state of all variables
#lwstate

# Calc unperturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)

unpert_lwup = lwdiagnostics['upwelling_longwave_flux_in_air']
unpert_lwdown = lwdiagnostics['downwelling_longwave_flux_in_air']


# prepare data to save to file
swup_data = np.reshape(unpert_swup.values,42)*0.965529262*0.501472063
swdown_data = np.reshape(unpert_swdown.values, 42)*0.965529262*0.501472063
lwup_data = np.reshape(unpert_lwup.values, 42) 
lwdown_data = np.reshape(unpert_lwdown.values, 42)
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)

# write the data out as json
with open("fluxes_330.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_data, 'sw_dn':swdown_data, 'lw_up':lwup_data, 'lw_dn':lwdown_data}, default=lambda x: list(x), indent=4))


