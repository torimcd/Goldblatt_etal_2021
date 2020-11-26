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
#print(swstate)

# Calc unperturbed fluxes - CO2 at 330 ppm default
swtendencies, swdiagnostics = swradiation(swstate)

unpert_swup = swdiagnostics['upwelling_shortwave_flux_in_air']
unpert_swdown = swdiagnostics['downwelling_shortwave_flux_in_air']

print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# ----------------------  Change CO2 concentration - 100ppm --------------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*0.3030303
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_100 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_100 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])


# prepare data to save to file
swup_100_data = np.reshape(swup_100,42)*0.965529262*0.501472063
swdown_100_data = np.reshape(swdown_100, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)




# ------------------ Change CO2 concentration - 500ppm ------------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*5
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_500 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_500 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_500_data = np.reshape(swup_500,42)*0.965529262*0.501472063
swdown_500_data = np.reshape(swdown_500, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





# --------------------------- Change CO2 concentration - 1000ppm ----------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*2
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_1000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_1000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_1000_data = np.reshape(swup_1000,42)*0.965529262*0.501472063
swdown_1000_data = np.reshape(swdown_1000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





# --------------------- Change CO2 concentration - 1500ppm -------------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.5
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_1500 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_1500 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_1500_data = np.reshape(swup_1500,42)*0.965529262*0.501472063
swdown_1500_data = np.reshape(swdown_1500, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





# ------------------- Change CO2 concentration - 2000ppm -----------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.333333333
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_2000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_2000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_2000_data = np.reshape(swup_2000,42)*0.965529262*0.501472063
swdown_2000_data = np.reshape(swdown_2000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





# ------------------------- Change CO2 concentration - 5000ppm -----------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*2.5
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_5000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_5000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_5000_data = np.reshape(swup_5000,42)*0.965529262*0.501472063
swdown_5000_data = np.reshape(swdown_5000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





#  -------------------- Change CO2 concentration - 10000ppm --------------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*2
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_10000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_10000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_10000_data = np.reshape(swup_10000,42)*0.965529262*0.501472063
swdown_10000_data = np.reshape(swdown_10000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





# ------------------- Change CO2 concentration - 20000ppm --------------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*2
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_20000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_20000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_20000_data = np.reshape(swup_20000,42)*0.965529262*0.501472063
swdown_20000_data = np.reshape(swdown_20000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)



# -------------------- Change CO2 concentration - 30000ppm -------------------------------------
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.5
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_30000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_30000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_30000_data = np.reshape(swup_30000,42)*0.965529262*0.501472063
swdown_30000_data = np.reshape(swdown_30000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)





# Change CO2 concentration - 35000ppm
swstate['mole_fraction_of_carbon_dioxide_in_air'].values = swstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.166666667
print(swstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
swtendencies, swdiagnostics = swradiation(swstate)
swup_35000 = np.array(swdiagnostics['upwelling_shortwave_flux_in_air'])
swdown_35000 = np.array(swdiagnostics['downwelling_shortwave_flux_in_air'])

# prepare data to save to file
swup_35000_data = np.reshape(swup_35000,42)*0.965529262*0.501472063
swdown_35000_data = np.reshape(swdown_35000, 42)*0.965529262*0.501472063
c5_pressure = np.reshape(swstate['air_pressure_on_interface_levels'].values, 42)




# ---------------------------------  Longwave  ---------------------------------------------------
lwradiation = climt.RRTMGLongwave()
surface = climt.SlabSurface()

# Create longwave state
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

unpert_lwup = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
unpert_lwdown = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# -------------------   Change CO2 concentration - 100ppm  -----------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*0.3030303
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_100 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_100 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_100_data = np.reshape(lwup_100, 42) 
lwdown_100_data = np.reshape(lwdown_100, 42)



# --------------------  Change CO2 concentration - 500ppm  --------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*5
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_500 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_500 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_500_data = np.reshape(lwup_500, 42) 
lwdown_500_data = np.reshape(lwdown_500, 42)



# --------------------- Change CO2 concentration - 1000ppm -----------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*2
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_1000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_1000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_1000_data = np.reshape(lwup_1000, 42) 
lwdown_1000_data = np.reshape(lwdown_1000, 42)






# ---------------------------- Change CO2 concentration - 1500ppm ----------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.5
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_1500 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_1500 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_1500_data = np.reshape(lwup_1500, 42) 
lwdown_1500_data = np.reshape(lwdown_1500, 42)




# ---------------------------- Change CO2 concentration - 2000ppm --------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.3333333333
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_2000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_2000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_2000_data = np.reshape(lwup_2000, 42) 
lwdown_2000_data = np.reshape(lwdown_2000, 42)





# ---------------------------- Change CO2 concentration - 5000ppm -------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*2.5
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_5000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_5000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_5000_data = np.reshape(lwup_5000, 42) 
lwdown_5000_data = np.reshape(lwdown_5000, 42)



 
# --------------------------  Change CO2 concentration - 10000ppm ----------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*2
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_10000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_10000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_10000_data = np.reshape(lwup_10000, 42) 
lwdown_10000_data = np.reshape(lwdown_10000, 42)




# --------------------------- Change CO2 concentration - 20000ppm -----------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*2
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_20000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_20000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_20000_data = np.reshape(lwup_20000, 42) 
lwdown_20000_data = np.reshape(lwdown_20000, 42)





# ------------------------ Change CO2 concentration - 30000ppm ---------------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.5
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_30000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_30000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_30000_data = np.reshape(lwup_30000, 42) 
lwdown_30000_data = np.reshape(lwdown_30000, 42)




# -------------------------- Change CO2 concentration - 35000ppm --------------------------
lwstate['mole_fraction_of_carbon_dioxide_in_air'].values = lwstate['mole_fraction_of_carbon_dioxide_in_air'].values*1.166666666667
print(lwstate['mole_fraction_of_carbon_dioxide_in_air'].values)

# Calc perturbed fluxes
lwtendencies, lwdiagnostics = lwradiation(lwstate)
lwup_35000 = np.array(lwdiagnostics['upwelling_longwave_flux_in_air'])
lwdown_35000 = np.array(lwdiagnostics['downwelling_longwave_flux_in_air'])

# prepare data to save to file
lwup_35000_data = np.reshape(lwup_35000, 42) 
lwdown_35000_data = np.reshape(lwdown_35000, 42)



# write the data out as json
with open("fluxes_100.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_100_data, 'sw_dn':swdown_100_data, 'lw_up':lwup_100_data, 'lw_dn':lwdown_100_data}, default=lambda x: list(x), indent=4))


with open("fluxes_500.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_500_data, 'sw_dn':swdown_500_data, 'lw_up':lwup_500_data, 'lw_dn':lwdown_500_data}, default=lambda x: list(x), indent=4))


with open("fluxes_1000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_1000_data, 'sw_dn':swdown_1000_data, 'lw_up':lwup_1000_data, 'lw_dn':lwdown_1000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_1500.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_1500_data, 'sw_dn':swdown_1500_data, 'lw_up':lwup_1500_data, 'lw_dn':lwdown_1500_data}, default=lambda x: list(x), indent=4))


with open("fluxes_2000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_2000_data, 'sw_dn':swdown_2000_data, 'lw_up':lwup_2000_data, 'lw_dn':lwdown_2000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_5000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_5000_data, 'sw_dn':swdown_5000_data, 'lw_up':lwup_5000_data, 'lw_dn':lwdown_5000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_10000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_10000_data, 'sw_dn':swdown_10000_data, 'lw_up':lwup_10000_data, 'lw_dn':lwdown_10000_data}, default=lambda x: list(x), indent=4))


with open("fluxes_20000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_20000_data, 'sw_dn':swdown_20000_data, 'lw_up':lwup_20000_data, 'lw_dn':lwdown_20000_data}, default=lambda x: list(x), indent=4))



with open("fluxes_30000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_30000_data, 'sw_dn':swdown_30000_data, 'lw_up':lwup_30000_data, 'lw_dn':lwdown_30000_data}, default=lambda x: list(x), indent=4))



with open("fluxes_35000.txt", 'w+') as datafile_id:
	datafile_id.write(json.dumps({'pressure':c5_pressure, 'sw_up':swup_35000_data, 'sw_dn':swdown_35000_data, 'lw_up':lwup_35000_data, 'lw_dn':lwdown_35000_data}, default=lambda x: list(x), indent=4))
