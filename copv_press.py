
# Regulated pressure-fed fluids system COPV pressure calculator

# This program aims to simulate the flow of pressurant gas through the dome regulators
# of a pressure-fed rocket fluids system in order to determine the initial amount of 
# pressure required to maintain pressure throughout the entire burn.

import numpy as np
import math

# Universal constants

R = 1.985875279009 # BTU / (lbmol * R)

# System definitions

copv_vol = 732 # cu. in.
ox_tank_vol = 924 # cu. in.
fuel_tank_vol = 924 # cu. in.

ox_dr_press = 570 # psi, ox side dome regulator pressure setting
fuel_dr_press = 630 # psi, fuel side dome regulator pressure setting

T_amb = 530 # R, ambient temperature

copv_press = 3500 # psi, initial COPV pressure

ox_press = 0 # psig, initial ox tank pressure
fuel_press = 0 # psig, initial fuel tank pressure

ox_vol = 100 # cu. in., initial gas volume in ox tank
fuel_vol = 100 # cu. in., initial gas volume in fuel tank

ox_vol_rate = -0.0376 * 12**3 # in^3/s, liquid oxidizer volume flow rate
fuel_vol_rate = -0.0291 * 12**3 # in^3/s, liquid fuel volume flow rate

# Liquid oxidizer properties

ox_name = "Liquid oxygen"

ox_M = 32 # lbm / lbmol, molar mass of oxidizer
T_ox_boil = 162.4 # R, oxidizer boiling point
T_ox = T_ox_boil # R, oxidizer initial temperature
ox_rho_l = 71.16 / 12**3 # lb/in^3
ox_heat_vap = 92 # BTU / lbm
print(ox_heat_vap * ox_M / (R*(1/170 - 1/T_ox_boil)))
ox_P_part = 14.7 / (np.exp(ox_heat_vap * ox_M / (R * (1/170 - 1/T_ox_boil))))
print(ox_P_part)

# Liquid fuel properties

fuel_name = "Kerosene"

fuel_M = 177 # lbm / lbmol, molar mass of fuel
T_fuel = T_amb # R, fuel initial temperature
T_fuel_boil = 770 # R, fuel boiling point
fuel_rho_l = 51.1 / 12**3 # lb/in^3

# Ullage gas properties

ul_name = "Nitrogen"

ul_M = 28 # lbm / lbmol, molar mass of pressurant gas
 
