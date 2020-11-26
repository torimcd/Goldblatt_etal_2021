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




# ----------------------  Change CO2 concentration - 100ppm --------------------------
state['mole_fraction_of_carbon_dioxide_in_air'].values = state['mole_fraction_of_carbon_dioxide_in_air'].values*0.3030303

print(state['mole_fraction_of_carbon_dioxide_in_air'])

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





# ------------------------- Change CO2 concentration - 35000ppm-------------------------------- 
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




# write the data out as json
with open("fluxes_100.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_100_data, 'sw_dn':swdown_100_data, 'lw_up':lwup_100_data, 'lw_dn':lwdown_100_data, 'net':net_flux_100_data}, default=lambda x: list(x), indent=4))

with open("fluxes_330.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_330_data, 'sw_dn':swdown_330_data, 'lw_up':lwup_330_data, 'lw_dn':lwdown_330_data, 'net':net_flux_330_data}, default=lambda x: list(x), indent=4))

with open("fluxes_500.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_500_data, 'sw_dn':swdown_500_data, 'lw_up':lwup_500_data, 'lw_dn':lwdown_500_data, 'net':net_flux_500_data}, default=lambda x: list(x), indent=4))


with open("fluxes_1000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_1000_data, 'sw_dn':swdown_1000_data, 'lw_up':lwup_1000_data, 'lw_dn':lwdown_1000_data, 'net':net_flux_1000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_1500.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_1500_data, 'sw_dn':swdown_1500_data, 'lw_up':lwup_1500_data, 'lw_dn':lwdown_1500_data, 'net':net_flux_1500_data}, default=lambda x: list(x), indent=4))


with open("fluxes_2000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_2000_data, 'sw_dn':swdown_2000_data, 'lw_up':lwup_2000_data, 'lw_dn':lwdown_2000_data, 'net':net_flux_2000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_5000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_5000_data, 'sw_dn':swdown_5000_data, 'lw_up':lwup_5000_data, 'lw_dn':lwdown_5000_data, 'net':net_flux_5000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_10000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_10000_data, 'sw_dn':swdown_10000_data, 'lw_up':lwup_10000_data, 'lw_dn':lwdown_10000_data, 'net':net_flux_10000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_20000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_20000_data, 'sw_dn':swdown_20000_data, 'lw_up':lwup_20000_data, 'lw_dn':lwdown_20000_data, 'net':net_flux_20000_data}, default=lambda x: list(x), indent=4))



with open("fluxes_30000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_30000_data, 'sw_dn':swdown_30000_data, 'lw_up':lwup_30000_data, 'lw_dn':lwdown_30000_data, 'net':net_flux_30000_data}, default=lambda x: list(x), indent=4))



with open("fluxes_35000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_35000_data, 'sw_dn':swdown_35000_data, 'lw_up':lwup_35000_data, 'lw_dn':lwdown_35000_data, 'net':net_flux_35000_data}, default=lambda x: list(x), indent=4))
